#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import os
import pickle
import shutil
import sqlite3
import struct
import zlib
from multiprocessing import Process, Pool, Queue
from tempfile import mkdtemp, mkstemp
from typing import Any, Callable, Generator, Iterable, List, Tuple, Union


def serialize(value: dict) -> bytes:
    return pickle.dumps(value, pickle.HIGHEST_PROTOCOL)


def traverse(src: dict, dst: dict):
    for k, v in src.items():
        if k in dst:
            if isinstance(v, dict):
                traverse(v, dst[k])
            elif isinstance(v, (list, tuple)):
                dst[k] += v
            elif isinstance(v, set):
                dst[k] |= v
            else:
                dst[k] = v
        else:
            dst[k] = v


class Bucket(object):
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.data = {}

    def __setitem__(self, key: Union[str, int], value: Any):
        self.data[key] = value

    def __iter__(self):
        with open(self.filepath, "rb") as fh:
            while True:
                try:
                    n_bytes, = struct.unpack("<L", fh.read(4))
                except struct.error:
                    break

                chunk = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                for k, v in chunk.items():
                    yield k, v

    @property
    def size(self) -> int:
        return os.path.getsize(self.filepath)

    def append(self, key: Union[str, int], value: Any):
        if key in self.data:
            self.data[key].append(value)
        else:
            self.data[key] = [value]

    def add(self, key: Union[str, int], value: Any):
        if key in self.data:
            self.data[key].add(value)
        else:
            self.data[key] = {value}

    def update(self, key: Union[str, int], value: dict):
        if key in self.data:
            traverse(value, self.data[key])
        else:
            self.data[key] = value

    def update_from_seq(self, key: Union[str, int], *args: Iterable):
        if key in self.data:
            d = self.data[key]
        else:
            d = self.data[key] = {} if len(args) > 1 else set()

        n = len(args) - 2
        for i, k in enumerate(args[:-1]):
            if k in d:
                d = d[k]
            else:
                d[k] = {} if i < n else set()
                d = d[k]

        d.add(args[-1])

    def sync(self):
        if self.data:
            s = zlib.compress(serialize(self.data))
            self.data = {}
            with open(self.filepath, "ab") as fh:
                fh.write(struct.pack("<L", len(s)) + s)

    def merge_item(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            data[k] = v

        return data

    def merge_list(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            if k in data:
                data[k] += v
            else:
                data[k] = v

        return data

    def merge_set(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            if k in data:
                data[k] |= v
            else:
                data[k] = v

        return data

    def merge_dict(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            if k in data:
                traverse(v, data[k])
            else:
                data[k] = v

        return data

    def merge(self, _type: Union[type, None]) -> dict:
        if _type is None:
            return self.merge_item()
        elif _type == list:
            return self.merge_list()
        elif _type == set:
            return self.merge_set()
        elif _type == dict:
            return self.merge_dict()
        else:
            raise ValueError(_type)


class Store(object):
    def __init__(self, filepath: str, keys: list=list(), processes: int=0,
                 tmpdir: Union[str, None]=None, cache_size: int=0,
                 dir_limit: int=1000):
        self.filepath = filepath
        self.keys = keys
        self.processes = processes

        # Root directory
        self.dir = None

        # Working directory (either root, or subdirectory)
        self._dir = None

        # Max number of files / directory
        self.dir_limit = dir_limit

        # Number of files in the current directory (self._dir)
        self.dir_count = 0

        # Multiprocessing variables
        self.workers = []               # pool of workers
        self.workers_keys = []          # keys assigned to each worker
        self.workers_items = []         # items to be feed to worker
        self.cache_size = cache_size    # max number of items in memory
        self.queues = []                # each worker has its own queue

        """
        Define which method from Bucket is to be used:
            0: __setitem__
            1: append
            2: add
            3: update
            4: update_from_seq
        """
        self.type = 0

        # List of Bucket instances, one per key
        self.buckets = []

        # Items currently in cache (when reading the Store)
        self.items = {}
        # File object (only when reading)
        self.fh = None

        # Define methods to use based on the mode (singe/multi-processing)
        if self.processes > 0:
            # Update cache size
            self.cache_size //= self.processes

            self._set_item = self._set_item_mp
            self.append = self._append_mp
            self.add = self._add_mp
            self.update = self._update_mp
            self.update_from_seq = self._update_from_seq_mp
            self.sync = self._sync_mp
            self.merge = self._merge_mp
            self._iter = self._iter_mp
        else:
            self._set_item = self._set_item_sp
            self.append = self._append_sp
            self.add = self._add_sp
            self.update = self._update_sp
            self.update_from_seq = self._update_from_seq_sp
            self.sync = self._sync_sp
            self.merge = self._merge_sp
            self._iter = self._iter_sp

        # Init based on mode (write/read)
        if self.keys:
            # Write mode

            if tmpdir is not None:
                os.makedirs(tmpdir, exist_ok=True)
            self._dir = self.dir = mkdtemp(dir=tmpdir)

            # Create one bucket per key
            self.buckets = [
                Bucket(self._mktemp())
                for _ in self.keys
            ]

            if self.processes > 0:
                # Multiprocessing mode
                step = int(len(self.keys) / self.processes + 1)

                for i in range(0, len(self.keys), step):
                    _keys = self.keys[i:i+step]
                    _buckets = self.buckets[i:i+step]

                    q = Queue(maxsize=1)
                    p = Process(
                        target=self._fill_buckets,
                        args=(_keys, _buckets, q)
                    )

                    self.workers.append(p)
                    self.workers_keys.append(self.keys[i])
                    self.workers_items.append([])
                    self.queues.append(q)

                    p.start()
        else:
            self.peek()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __setitem__(self, key, value):
        self._set_item(key, value)

    def __getitem__(self, key: Union[str, int]):
        if key in self.items:
            return self.items[key]

        i = bisect.bisect_right(self.keys, key)
        if not i:
            raise KeyError(key)

        try:
            offset = self.offsets[i-1]
        except IndexError:
            raise KeyError(key)

        if self.load_chunk(offset):
            return self.items[key]
        else:
            raise KeyError(key)

    def __iter__(self) -> Generator[Tuple, None, None]:
        return self._iter()

    @property
    def size(self) -> int:
        return sum([bucket.size for bucket in self.buckets])

    def close(self):
        if self.dir:
            shutil.rmtree(self.dir)
            self.dir = None

    def load_chunk(self, offset) -> bool:
        if self.offset == offset:
            return False
        elif self.fh is None:
            self.fh = open(self.filepath, "rb")

        self.fh.seek(offset)
        self.offset = offset

        n_bytes, = struct.unpack("<L", self.fh.read(4))
        self.items = pickle.loads(zlib.decompress(self.fh.read(n_bytes)))
        return True

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def peek(self) -> bool:
        try:
            fh = open(self.filepath, "rb")
        except FileNotFoundError:
            return False

        try:
            offset, = struct.unpack("<Q", fh.read(8))
        except struct.error:
            return False
        else:
            fh.seek(offset)
            try:
                keys, offsets = pickle.load(fh)
            except pickle.UnpicklingError:
                return False
            else:
                if not self.keys:
                    self.keys = keys
                if not self.offsets:
                    self.offsets = offsets
                return True
        finally:
            fh.close()

    def get_bucket(self, key: Union[str, int]) -> Bucket:
        i = bisect.bisect_right(self.keys, key)
        if i:
            return self.buckets[i-1]
        else:
            raise KeyError(key)

    def _assign_to_worker(self, key: Union[str, int], value: Any):
        i = bisect.bisect_right(self.workers_keys, key)
        if i:
            i -= 1
            self.workers_items[i].append((key, value))
            if len(self.workers_items[i]) == self.cache_size:
                self.queues[i].put((self.type, self.workers_items[i]))
                self.workers_items[i] = []
        else:
            raise KeyError(key)

    def get_type(self) -> Union[type, None]:
        if not self.type:
            return None
        elif self.type == 1:
            return list
        elif self.type == 2:
            return set
        elif self.type in (3, 4):
            return dict
        else:
            raise ValueError(self.type)

    def _set_item(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _set_item_sp(self, key: Union[str, int], value: Any):
        self.get_bucket(key)[key] = value
        self.type = 0

    def _set_item_mp(self, key: Union[str, int], value: Any):
        self._assign_to_worker(key, value)
        self.type = 0

    def append(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _append_sp(self, key: Union[str, int], value: Any):
        self.get_bucket(key).append(key, value)
        self.type = 1

    def _append_mp(self, key: Union[str, int], value: Any):
        self._assign_to_worker(key, value)
        self.type = 1

    def add(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _add_sp(self, key: Union[str, int], value: Any):
        self.get_bucket(key).add(key, value)
        self.type = 2

    def _add_mp(self, key: Union[str, int], value: Any):
        self._assign_to_worker(key, value)
        self.type = 2

    def update(self, key: Union[str, int], value: dict):
        raise NotImplementedError

    def _update_sp(self, key: Union[str, int], value: dict):
        self.get_bucket(key).update(key, value)
        self.type = 3

    def _update_mp(self, key: Union[str, int], value: dict):
        self._assign_to_worker(key, value)
        self.type = 3

    def update_from_seq(self, key: Union[str, int], *args: Iterable):
        raise NotImplementedError

    def _update_from_seq_sp(self, key: Union[str, int], *args: Iterable):
        self.get_bucket(key).update_from_seq(key, args)
        self.type = 4

    def _update_from_seq_mp(self, key: Union[str, int], *args: Iterable):
        self._assign_to_worker(key, value)
        self.type = 4

    def sync(self):
        raise NotImplementedError

    def _sync_sp(self):
        for bucket in self.buckets:
            bucket.sync()

    def _sync_mp(self):
        for i, chunk in enumerate(self.workers_items):
            self.queues[i].put((self.type, chunk))
            self.workers_items[i] = []

    def merge(self, func: Callable=None):
        raise NotImplementedError

    def _merge_sp(self, func: Callable=None):
        self.sync()
        pos = 0
        self.offsets = []

        _type = self.get_type()
        with open(self.filepath, "wb") as fh:
            # Header (empty for now)
            pos += fh.write(struct.pack("<Q", 0))

            # Body
            for bucket in self.buckets:
                items = bucket.merge(_type)

                if func is not None:
                    self._dapply(items, func)

                self.offsets.append(pos)
                chunk = zlib.compress(serialize(items))
                pos += fh.write(struct.pack("<L", len(chunk)) + chunk)

            # Footer
            pickle.dump((self.keys, self.offsets), fh)

            # Header
            fh.seek(0)
            fh.write(struct.pack("<Q", pos))

    def _merge_mp(self, func: Callable=None):
        self.sync()

        # Send signal to workers that we're done
        for q in self.queues:
            q.put(None)

        # Wait for workers to complete
        for w in self.workers:
            w.join()

        pos = 0
        self.offsets = []

        _type = self.get_type()
        with open(self.filepath, "wb") as fh, Pool(self.processes) as pool:
            # Header (empty for now)
            pos += fh.write(struct.pack("<Q", 0))

            iterable = [(bucket, _type, func) for bucket in self.buckets]
            for chunk in pool.imap(self._merge_bucket, iterable):
                self.offsets.append(pos)
                pos += fh.write(struct.pack("<L", len(chunk)) + chunk)

            # Footer
            pickle.dump((self.keys, self.offsets), fh)

            # Header
            fh.seek(0)
            fh.write(struct.pack("<Q", pos))

    def _iter(self) -> Generator[Tuple, None, None]:
        raise NotImplementedError

    def _iter_sp(self) -> Generator[Tuple, None, None]:
        for offset in self.offsets:
            self.load_chunk(offset)
            for key in sorted(self.items):
                yield key, self.items[key]

    def _iter_mp(self) -> Generator[Tuple, None, None]:
        with Pool(self.processes) as pool:
            iterable = [(self.filepath, offset) for offset in self.offsets]
            for items in pool.imap(self._load_chunk, iterable):
                for key, value in items:
                    yield key, value

    @staticmethod
    def _load_chunk(filepath: str, offset: str) -> list:
        with open(filepath, "rb") as fh:
            fh.seek(offset)
            n_bytes, = struct.unpack("<L", fh.read(4))
            items = pickle.loads(zlib.decompress(fh.read(n_bytes)))

        return [(key, items[key]) for key in sorted(items)]

    def _mktemp(self) -> str:
        if self.dir_count + 1 == self.dir_limit:
            # Too many files in directory: create a subdirectory
            self._dir = mkdtemp(dir=self._dir)
            self.dir_count = 0

        fd, filepath = mkstemp(dir=self._dir)
        os.close(fd)
        self.dir_count += 1
        return filepath

    @staticmethod
    def _fill_buckets(keys: List[str], buckets: List[Bucket], queue: Queue):
        while True:
            task = queue.get()
            if task is None:
                break

            _type, items = task
            for key, value in items:
                i = bisect.bisect_right(keys, key)
                if i:
                    bucket = buckets[i-1]

                    if not _type:
                        bucket[key] = value
                    elif _type == 1:
                        bucket.append(key, value)
                    elif _type == 2:
                        bucket.add(key, value)
                    elif _type == 3:
                        bucket.update(key, value)
                    elif _type == 4:
                        bucket.update_from_seq(key, value)
                    else:
                        raise ValueError(_type)
                else:
                    raise KeyError(key)

            for bucket in buckets:
                bucket.sync()

    @staticmethod
    def _dapply(data: dict, func: Callable):
        for k, v in data.items():
            data[k] = func(v)

    @staticmethod
    def _merge_bucket(args: Tuple[Bucket, Union[type, None],
                      Union[Callable, None]]) -> bytes:
        bucket, _type, func = args
        items = bucket.merge(_type)

        if func is not None:
            Store._dapply(items, func)

        return zlib.compress(serialize(items))


class KVdb(object):
    def __init__(self, filepath: str=None, cache_size: int=0):
        if filepath:
            self.filepath = filepath
            self.temporary = False
        else:
            fd, self.filepath = mkstemp()
            os.close(fd)
            os.remove(self.filepath)
            self.temporary = True

        self.con = sqlite3.connect(self.filepath)
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS data (
                id TEXT PRIMARY KEY NOT NULL,
                val TEXT NOT NULL
            )
            """
        )
        self.cache_size = cache_size
        self.cache_items = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __setitem__(self, key: str, value: Any):
        if key in self.cache_items:
            self.cache_items[key] = value
        elif self.cache_size:
            if len(self.cache_items) == self.cache_size:
                self.sync()
            self.cache_items[key] = value
        else:
            self.con.execute(
                "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)",
                (key, serialize(value)
            )
            self.con.commit()

    def __getitem__(self, key: str) -> dict:
        if key in self.cache_items:
            return self.cache_items[key]
        else:
            cur = self.con.execute("SELECT val FROM data WHERE id=?", (key,))
            row = cur.fetchone()
            if row is None:
                raise KeyError(key)

            value = pickle.loads(row[0])
            if self.cache_size:
                if len(self.cache_items) == self.cache_size:
                    self.sync()
                self.cache_items[key] = value

            return value

    def __iter__(self) -> Generator:
        self.sync()

        # Disable cache because each key/value pair is accessed once
        cache_size = self.cache_size
        self.cache_size = 0

        keys = [row[0] for row in self.con.execute("SELECT id FROM data")]
        for key in keys:
            yield key, self[key]

        # Restore user-defined cache size
        self.cache_size = cache_size

    def sync(self):
        if self.cache_items:
            self.con.executemany(
                "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)",
                (
                    (key, serialize(value)
                    for key, value in self.cache_items.items()
                )
            )
            self.con.commit()
            self.cache_items = {}

    def close(self):
        if self.con is not None:
            self.sync()
            self.con.close()
            self.con = None

        if self.filepath and self.temporary:
            os.remove(self.filepath)
            self.filepath = None

    def getsize(self) -> int:
        try:
            return os.path.getsize(self.filepath)
        except (FileNotFoundError, TypeError):
            return 0
