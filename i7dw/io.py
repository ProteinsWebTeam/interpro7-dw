#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import logging
import os
import pickle
import shutil
import sqlite3
import struct
import zlib
from multiprocessing import Process, Pool, Queue
from tempfile import mkdtemp, mkstemp
from typing import Any, Callable, Generator, Iterable, Tuple, Union


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


class Shelf(object):
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

    def flush(self):
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

    def merge(self, _type: type) -> Tuple[dict, int]:
        if _type == list:
            items = self.merge_list()
        elif _type == set:
            items = self.merge_set()
        elif _type == dict:
            items = self.merge_dict()
        else:
            items = self.merge_item()

        return items, os.path.getsize(self.filepath)


class Aisle(object):
    def __init__(self, keys: list, dir: str=None, dir_limit: int=1000):
        self.keys = keys

        # If not None: expects existing
        self.dir = mkdtemp(dir=dir)
        self._dir = self.dir

        self.dir_limit = dir_limit
        self.dir_count = 0
        self.shelves = [self.create_shelf() for _ in self.keys]

        """
        Type of values:
            * None: overwrite existing items
            * list: append to list
            * set: add to set
            * dict: update dict
                (if a sub-key of the dict already exists,
                the value will be updated accordingly to its type)
        """
        self.type = None

    def create_shelf(self) -> Shelf:
        if self.dir_count + 1 == self.dir_limit:
            # Too many files in directory: create a subdirectory
            self._dir = mkdtemp(dir=self._dir)
            self.dir_count = 0

        fd, filepath = mkstemp(dir=self._dir)
        os.close(fd)
        self.dir_count += 1
        return Shelf(filepath)

    def __setitem__(self, key: Union[str, int], value: Any):
        shelf = self.get_shelf(key)
        shelf[key] = value
        self.type = None

    def append(self, key: Union[str, int], value: Any):
        shelf = self.get_shelf(key)
        shelf.append(key, value)
        self.type = list

    def add(self, key: Union[str, int], value: Any):
        shelf = self.get_shelf(key)
        shelf.add(key, value)
        self.type = set

    def update(self, key: Union[str, int], value: dict):
        shelf = self.get_shelf(key)
        shelf.update(key, value)
        self.type = dict

    def update_from_seq(self, key: Union[str, int], args: Iterable):
        shelf = self.get_shelf(key)
        shelf.update_from_seq(key, *args)
        self.type = dict

    def get_shelf(self, key: Union[str, int]) -> Shelf:
        i = bisect.bisect_right(self.keys, key)
        if i :
            return self.shelves[i-1]
        else:
            raise KeyError(key)

    def flush(self):
        for shelf in self.shelves:
            shelf.flush()

    def merge(self) -> Generator[Tuple[dict, int], None, None]:
        i = 0
        d = os.path.basename(self.dir)
        for shelf in self.shelves:
            yield shelf.merge(self.type)
            i += 1
            logging.info("{} merged {}".format(d, i))


class Store2(object):
    def __init__(self, filepath: str, keys: list=list(), processes: int=0,
                dir: str=None):
        self.filepath = filepath
        self.keys = keys
        self.processes = processes
        self.dir = dir
        if self.dir:
            os.makedirs(self.dir, exist_ok=True)

        self.queue_in = None
        self.queues_out = []
        self.chunk = []
        self.aisles = []
        if self.keys:
            if self.processes > 0:
                self.queue_in = Queue(self.processes)
                for _ in range(self.processes):
                    queue_out = Queue(1)
                    p = Process(target=self._create_aisle,
                                args=(self.keys, self.dir, self.queue_in,
                                      queue_out)
                    )
                    p.start()
                    self.queues_out.append(queue_out)
                    self.aisles.append(p)

                self._set_item = self._set_item_mp
                self.append = self._append_mp
                self.add = self._add_mp
                self.update = self._update_mp
                self.update_from_seq = self._update_from_seq_mp
                self.flush = self._flush_mp
                self.merge = self._merge_mp
            else:
                self._set_item = self._set_item_sp
                self.aisles.append(Aisle(self.keys, self.dir))

                self._set_item = self._set_item_sp
                self.append = self._append_sp
                self.add = self._add_sp
                self.update = self._update_sp
                self.update_from_seq = self._update_from_seq_sp
                self.flush = self._flush_sp
                self.merge = self._merge_sp

        self.offsets = []

        """
        Used only for multiprocessing
        Define which method from Aisle is to be used:
            0: __setitem__
            1: append
            2: add
            3: update
            4: update_from_seq
        """
        self.type = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __setitem__(self, key: Union[str, int], value: Any):
        self._set_item(key, value)

    def close(self):
        pass

    def peek(self) -> bool:
        pass

    def _set_item(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _set_item_sp(self, key: Union[str, int], value: Any):
        self.aisles[0][key] = value

    def _set_item_mp(self, key: Union[str, int], value: Any):
        self.chunk.append((key, value))
        self.type = 0

    def append(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _append_sp(self, key: Union[str, int], value: Any):
        self.aisles[0].append(key, value)

    def _append_mp(self, key: Union[str, int], value: Any):
        self.chunk.append((key, value))
        self.type = 1

    def add(self, key: Union[str, int], value: Any):
        raise NotImplementedError

    def _add_sp(self, key: Union[str, int], value: Any):
        self.aisles[0].add(key, value)

    def _add_mp(self, key: Union[str, int], value: Any):
        self.chunk.append((key, value))
        self.type = 2

    def update(self, key: Union[str, int], value: dict):
        raise NotImplementedError

    def _update_sp(self, key: Union[str, int], value: dict):
        self.aisles[0].update(key, value)

    def _update_mp(self, key: Union[str, int], value: dict):
        self.chunk.append((key, value))
        self.type = 3

    def update_from_seq(self, key: Union[str, int], *args: Iterable):
        raise NotImplementedError

    def _update_from_seq_sp(self, key: Union[str, int], *args: Iterable):
        self.aisles[0].update_from_seq(key, args)

    def _update_from_seq_mp(self, key: Union[str, int], *args: Iterable):
        self.chunk.append((key, args))
        self.type = 4

    def flush(self):
        raise NotImplementedError

    def _flush_sp(self):
        self.aisles[0].flush()

    def _flush_mp(self):
        self.queue_in.put((self.type, self.chunk))
        self.chunk = []

    def merge(self, func: Callable=None) -> int:
        raise NotImplementedError

    def _merge_sp(self, func: Callable=None) -> int:
        self.flush()
        pos = 0
        size = 0
        self.offsets = []
        with open(self.filepath, "wb") as fh:
            pos += fh.write(struct.pack("<Q", 0))

            for items, _size in self.aisles[0].merge():
                size += _size

                if func is not None:
                    self.post(items, func)

                self.offsets.append(pos)
                chunk = zlib.compress(serialize(items))
                pos += fh.write(struct.pack("<L", len(chunk)) + chunk)

            fh.seek(0)
            fh.write(struct.pack("<Q", pos))

        return size

    def _merge_mp(self, func: Callable=None) -> int:
        self.flush()
        pos = 0
        size = 0
        self.offsets = []

        for _ in self.aisles:
            self.queue_in.put(None)

        with open(self.filepath, "wb") as fh:
            pos += fh.write(struct.pack("<Q", 0))

            i = 0
            for _ in self.keys:
                i += 1
                logging.info("store: {}: first aisle".format(i))
                items, _size = self.queues_out[0].get()
                size += _size

                for q in self.queues_out[1:]:
                    logging.info("store: {}: get from aisle".format(i))
                    _items, _size = q.get()
                    logging.info("store: {}: traverse".format(i))
                    size += _size
                    traverse(_items, items)

                logging.info("store: {}: apply func".format(i))
                if func is not None:
                    self.post(items, func)

                self.offsets.append(pos)
                logging.info("store: {}: compress".format(i))
                chunk = zlib.compress(serialize(items))
                logging.info("store: {}: write".format(i))
                pos += fh.write(struct.pack("<L", len(chunk)) + chunk)

            fh.seek(0)
            fh.write(struct.pack("<Q", pos))
        return size

    @staticmethod
    def post(data: dict, func: Callable):
        for k, v in data.items():
            data[k] = func(v)

    @staticmethod
    def _create_aisle(keys: list, dir: Union[str, None], queue_in: Queue,
                      queue_out: Queue):
        aisle = Aisle(keys, dir)

        types = {
            0: aisle.__setitem__,
            1: aisle.append,
            2: aisle.add,
            3: aisle.update,
            4: aisle.update_from_seq
        }

        while True:
            chunk = queue_in.get()
            if chunk is None:
                break

            _type, items = chunk
            func = types[_type]
            for key, value in items:
                func(key, value)

            aisle.flush()

        for items, size in aisle.merge():
            queue_out.put((items, size))


class Bucket(object):
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.data = {}

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

    def __setitem__(self, key: Union[str, int], value: Any):
        self.data[key] = value

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

    def update_from_seq(self, key: str, *args: Iterable):
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

    def update(self, key: Union[str, int], value: dict):
        if key in self.data:
            self.traverse(value, self.data[key])
        else:
            self.data[key] = value

    @staticmethod
    def traverse(src: dict, dst: dict):
        for k, v in src.items():
            if k in dst:
                if isinstance(v, dict):
                    Bucket.traverse(v, dst[k])
                elif isinstance(v, (list, tuple)):
                    dst[k] += v
                elif isinstance(v, set):
                    dst[k] |= v
                else:
                    dst[k] = v
            else:
                dst[k] = v

    def flush(self):
        if self.data:
            s = zlib.compress(pickle.dumps(self.data,
                                           pickle.HIGHEST_PROTOCOL))

            with open(self.filepath, "ab") as fh:
                fh.write(struct.pack("<L", len(s)) + s)

            self.data = {}

    def merge_dict(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            if k in data:
                self.traverse(v, data[k])
            else:
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

    def load(self) -> dict:
        data = self.data
        self.data = {}
        for k, v in self:
            data[k] = v

        return data

    def merge(self, _type: type):
        if _type == dict:
            return self.merge_dict()
        elif _type == list:
            return self.merge_list()
        elif _type == set:
            return self.merge_set()
        else:
            return self.load()


class Store(object):
    def __init__(self, filepath: str, keys: list=list(), processes: int=0,
                 tmpdir=None):
        self.filepath = filepath

        # To find chunk in filepath
        self.keys = keys
        self.offsets = []

        self.processes = processes

        if tmpdir is not None:
            os.makedirs(tmpdir, exist_ok=True)

        # Type of values stored (None, default: overwrite any existing value)
        self.type = None

        # Variables used when reading the file
        self.fh = None
        self.items = {}
        self.offset = None

        self.dir_limit = 1000
        self.dir_count = 0

        if self.keys:
            self.dir = mkdtemp(dir=tmpdir)
            self._dir = self.dir

            # Buckets when creating the file
            self.buckets = [self.create_bucket() for _ in self.keys]
        else:
            self._dir = self.dir = None
            self.buckets = []

        # Processes, used in iter()
        self.workers = None

        self.peek()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, key: Union[str, int]):
        if key in self.items:
            return self.items[key]

        i = bisect.bisect_right(self.keys, key)
        if i == 0:
            raise KeyError(key)

        try:
            offset = self.offsets[i-1]
        except IndexError:
            raise KeyError(key)

        if self.load_chunk(offset):
            return self.items[key]
        else:
            raise KeyError(key)

    def __setitem__(self, key: Union[str, int], value: Any):
        b = self.get_bucket(key)
        b[key] = value
        self.type = None

    def __iter__(self) -> Generator:
        for offset in self.offsets:
            self.load_chunk(offset)
            for key in sorted(self.items):
                yield key, self.items[key]

    def peek(self) -> bool:
        try:
            fh = open(self.filepath, "rb")
        except FileNotFoundError:
            return False

        try:
            offset, = struct.unpack('<Q', fh.read(8))
        except struct.error:
            return False
        else:
            fh.seek(offset)
            try:
                keys, offsets = pickle.loads(fh.read())
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

    def iter(self, processes: int=1) -> Generator:
        if processes > 1:
            return self._iter(processes-1)
        else:
            return self

    def _iter(self, processes: int=1) -> Generator:
        queue_in = Queue()
        queue_out = Queue(maxsize=processes)

        self.workers = [
            Process(target=self._load_chunk,
                    args=(self.filepath, queue_in, queue_out))
            for _ in range(processes)
        ]

        for w in self.workers:
            w.start()

        for offset in self.offsets:
            queue_in.put(offset)

        for _ in self.workers:
            queue_in.put(None)

        running_workers = len(self.workers)
        while True:
            items = queue_out.get()
            if items is None:
                running_workers -= 1
                if running_workers:
                    continue
                else:
                    break

            for key, value in items:
                yield key, value

        for w in self.workers:
            w.join()

    @staticmethod
    def _load_chunk(filepath: str, queue_in: Queue, queue_out: Queue):
        with open(filepath, "rb") as fh:
            while True:
                offset = queue_in.get()
                if offset is None:
                    break

                fh.seek(offset)
                n_bytes, = struct.unpack("<L", fh.read(4))
                items = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                queue_out.put([(key, items[key]) for key in sorted(items)])

        queue_out.put(None)

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

    def create_bucket(self):
        if self.dir_count + 1 == self.dir_limit:
            # Too many files in directory: create a subdirectory
            self._dir = mkdtemp(dir=self._dir)
            self.dir_count = 0

        fd, filepath = mkstemp(dir=self._dir)
        os.close(fd)
        self.dir_count += 1
        return Bucket(filepath)

    def get_bucket(self, key: Union[str, int]) -> Bucket:
        i = bisect.bisect_right(self.keys, key)
        if i:
            return self.buckets[i-1]
        else:
            raise ValueError("invalid key '{}'".format(key))

    def append(self, key, value):
        b = self.get_bucket(key)
        b.append(key, value)
        self.type = list

    def update(self, key: str, value: dict):
        b = self.get_bucket(key)
        b.update(key, value)
        self.type = dict

    def update_from_seq(self, key: str, *args: Iterable):
        b = self.get_bucket(key)
        b.update_from_seq(key, *args)
        self.type = dict

    def flush(self):
        for b in self.buckets:
            b.flush()

    @staticmethod
    def post(data: dict, func: Callable):
        for k, v in data.items():
            data[k] = func(v)

    @staticmethod
    def load_bucket(args):
        bucket, _type, func = args
        chunk = bucket.merge(_type)
        if func:
            Store.post(chunk, func)
        return chunk

    def close(self):
        self.items = {}
        self.offset = None

        if self.fh is not None:
            self.fh.close()
            self.fh = None

        if self.dir is not None:
            shutil.rmtree(self.dir)
            self.dir = None

        if self.workers is not None:
            for w in self.workers:
                if w.is_alive():
                    w.terminate()

            self.workers = None

    @staticmethod
    def merge_bucket(args: Tuple[Bucket, type, Callable]):
        b, _type, func = args
        data = b.merge(_type)
        os.remove(b.filepath)

        if func is not None:
            Store.post(data, func)

        return zlib.compress(pickle.dumps(data, pickle.HIGHEST_PROTOCOL))

    def merge_buckets(self, func: Callable=None,
                      processes: int=1) -> Generator:
        if processes > 1:
            with Pool(processes-1) as pool:
                iterable = [(b, self.type, func) for b in self.buckets]
                for chunk in pool.imap(self.merge_bucket, iterable):
                    yield chunk

        else:
            for b in self.buckets:
                data = b.merge(self.type)
                os.remove(b.filepath)

                if func is not None:
                    self.post(data, func)

                yield zlib.compress(pickle.dumps(data,
                                                 pickle.HIGHEST_PROTOCOL))

    def merge(self, func: Callable=None, processes: int=1) -> int:
        size = sum([os.path.getsize(b.filepath) for b in self.buckets])
        if self.buckets:
            self.flush()
            pos = 0
            self.offsets = []
            with open(self.filepath, "wb") as fh:
                pos += fh.write(struct.pack('<Q', 0))

                for chunk in self.merge_buckets(func, processes):
                    self.offsets.append(pos)
                    pos += fh.write(struct.pack("<L", len(chunk)) + chunk)

                fh.write(pickle.dumps((self.keys, self.offsets),
                                      pickle.HIGHEST_PROTOCOL))

                fh.seek(0)
                fh.write(struct.pack('<Q', pos))

            self.buckets = []

        return size

    def getsize(self) -> int:
        return sum([os.path.getsize(b.filepath) for b in self.buckets])

    def save(self):
        self.flush()
        with open(self.filepath + ".save", "wb") as fh:
            pickle.dump((self.dir,
                         self.keys,
                         self.type,
                         [b.filepath for b in self.buckets]), fh)

        self.dir = None

    def reload(self):
        filepath = self.filepath + ".save"
        with open(filepath, "rb") as fh:
            self.dir, self.keys, self.type, files = pickle.load(fh)

        self.buckets = [Bucket(f) for f in files]
        os.remove(filepath)


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
                (key, self.serialize(value))
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
        keys = [row[0] for row in self.con.execute("SELECT id FROM data")]
        for key in keys:
            yield key, self[key]

    def sync(self):
        if self.cache_items:
            self.con.executemany(
                "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)",
                (
                    (key, self.serialize(value))
                    for key, value in self.cache_items.items()
                )
            )
            self.con.commit()
            self.cache_items = {}

    @staticmethod
    def serialize(value: dict) -> bytes:
        return pickle.dumps(value, pickle.HIGHEST_PROTOCOL)

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
