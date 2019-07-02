# -*- coding: utf-8 -*-

import bisect
import os
import json
import pickle
import shutil
import sqlite3
import struct
import zlib
from multiprocessing import Pool, Queue
from tempfile import mkdtemp, mkstemp
from typing import Any, Callable, Generator, Iterable, Optional, Tuple, Union


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
    def __init__(self, filepath: Optional[str]=None, keys: list=list(),
                 tmpdir: Optional[str]=None, dir_limit: int=1000):
        self.filepath = filepath
        self.keys = keys

        # True if filepath is None (delete on close)
        self.temporary = False

        # Root directory
        self.dir = None

        # Working directory (either root, or subdirectory)
        self._dir = None

        # Max number of files / directory
        self.dir_limit = dir_limit

        # Number of files in the current directory (self._dir)
        self.dir_count = 0

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

        # Offsets of buckets (set in peek())
        self.offsets = []

        # Offset of currently loaded bucket
        self.offset = None

        # Init based on mode (write/read)
        if self.keys:
            # Write mode

            if tmpdir:
                os.makedirs(tmpdir, exist_ok=True)
            self._dir = self.dir = mkdtemp(dir=tmpdir)

            if self.filepath is None:
                self.temporary = True
                fd, self.filepath = mkstemp(dir=self.dir)
                os.close(fd)

            # Create one bucket per key
            self.buckets = [
                Bucket(self._mktemp())
                for _ in self.keys
            ]
        else:
            self.peek()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __setitem__(self, key: Union[str, int], value: Any):
        self.get_bucket(key)[key] = value
        self.type = 0

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
        for offset in self.offsets:
            self.load_chunk(offset)
            for key in sorted(self.items):
                yield key, self.items[key]

    @property
    def size(self) -> int:
        size = sum([bucket.size for bucket in self.buckets])
        if self.temporary:
            try:
                size += os.path.getsize(self.filepath)
            except FileNotFoundError:
                pass
        return size

    def close(self):
        self.items = {}

        if self.dir:
            shutil.rmtree(self.dir)
            self.dir = None

        if self.fh is not None:
            self.fh.close()
            self.fh = None

        if self.temporary:
            try:
                os.remove(self.filepath)
            except FileNotFoundError:
                pass

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

    def append(self, key: Union[str, int], value: Any):
        self.get_bucket(key).append(key, value)
        self.type = 1

    def add(self, key: Union[str, int], value: Any):
        self.get_bucket(key).add(key, value)
        self.type = 2

    def update(self, key: Union[str, int], value: dict):
        self.get_bucket(key).update(key, value)
        self.type = 3

    def update_from_seq(self, key: Union[str, int], *args: Iterable):
        self.get_bucket(key).update_from_seq(key, *args)
        self.type = 4

    def sync(self):
        for bucket in self.buckets:
            bucket.sync()

    def merge(self, func: Callable=None, processes: int=1) -> int:
        # Sync remaining items
        self.sync()

        size = self.size
        if processes > 1:
            self._merge_mp(processes-1, func)
        else:
            self._merge_sp(func)
        return max(size, self.size)

    def _merge_sp(self, func: Callable=None):
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

    def _merge_mp(self, processes:int, func: Callable=None):
        pos = 0
        self.offsets = []

        _type = self.get_type()
        with open(self.filepath, "wb") as fh, Pool(processes) as pool:
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

    @staticmethod
    def _load_chunk(filepath: str, input: Queue, output: Queue):
        for i, offset in iter(input.get, None):
            with open(filepath, "rb") as fh:
                fh.seek(offset)
                n_bytes, = struct.unpack("<L", fh.read(4))
                items = pickle.loads(zlib.decompress(fh.read(n_bytes)))

            output.put((i, [(key, items[key]) for key in sorted(items)]))

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
    def __init__(self, filepath: Optional[str]=None, dir: Optional[str]=None,
                 writeback: bool=False, insertonly: bool=False):
        if filepath:
            self.filepath = filepath
            self.temporary = False
        else:
            if dir:
                os.makedirs(dir, exist_ok=True)

            fd, self.filepath = mkstemp(dir=dir)
            os.close(fd)
            os.remove(self.filepath)
            self.temporary = True

        self.writeback = writeback
        self.insertonly = insertonly

        self.con = sqlite3.connect(self.filepath)
        if self.insertonly:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS data (
                    id TEXT NOT NULL,
                    val TEXT NOT NULL
                )
                """
            )
            self.stmt = "INSERT INTO data (id, val) VALUES (?, ?)"
        else:
            self.con.execute(
                """
                CREATE TABLE IF NOT EXISTS data (
                    id TEXT PRIMARY KEY NOT NULL,
                    val TEXT NOT NULL
                )
                """
            )
            self.stmt = "INSERT OR REPLACE INTO data (id, val) VALUES (?, ?)"
        self.cache = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if self.temporary:
            try:
                os.remove(self.filepath)
            except FileNotFoundError:
                pass

    def __del__(self):
        self.close()
        if self.temporary:
            try:
                os.remove(self.filepath)
            except FileNotFoundError:
                pass

    def __setitem__(self, key: str, value: Any):
        if self.writeback:
            self.cache[key] = value
        else:
            self.con.execute(self.stmt, (key, serialize(value)))
            self.con.commit()

    def __getitem__(self, key: str) -> Any:
        try:
            value = self.cache[key]
        except KeyError:
            row = self.con.execute(
                "SELECT val FROM data WHERE id=?", (key,)
            ).fetchone()
            if row:
                value = pickle.loads(row[0])
                if self.writeback:
                    self.cache[key] = value
            else:
                raise KeyError(key)

        return value

    def __iter__(self) -> Generator:
        self.close()
        with sqlite3.connect(self.filepath) as con:
            for row in con.execute("SELECT id, val FROM data ORDER BY id"):
                yield row[0], pickle.loads(row[1])

    def sync(self):
        if not self.cache:
            return

        self.con.executemany(
            self.stmt,
            ((key, serialize(value)) for key, value in self.cache.items())
        )
        self.con.commit()
        self.cache = {}

    def index():
        if self.insertonly:
            self.con.execute("CREATE UNIQUE INDEX idx_data ON data (id)")

    def close(self):
        if self.con is None:
            return

        self.sync()
        self.index()
        self.con.close()
        self.con = None

    @property
    def size(self) -> int:
        return os.path.getsize(self.filepath)


class TempFile(object):
    def __init__(self, dir: Optional[str]=None):
        if dir is not None:
            os.makedirs(dir, exist_ok=True)

        fd, self.path = mkstemp(dir=dir)
        os.close(fd)

        self.fh = open(self.path, "wb")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        self.remove()

    def __del__(self):
        self.close()
        self.remove()

    def __iter__(self):
        self.close()
        with open(self.path, "rb") as fh:
            while True:
                try:
                    obj = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield obj

    @property
    def size(self) -> int:
        try:
            return os.path.getsize(self.path)
        except FileNotFoundError:
            return 0

    def write(self, obj):
        pickle.dump(obj, self.fh)

    def close(self):
        self.fh.close()

    def remove(self):
        try:
            os.remove(self.path)
        except FileNotFoundError:
            pass


class JsonFileOrganizer(object):
    def __init__(self, root: str, items_per_file: int=10000,
                 files_per_dir: int=1000, func: Callable=lambda x: x,
                 indent: Optional[int]=None):
        os.makedirs(root, exist_ok=True)
        os.chmod(root, 0o775)
        self.root = root
        self.dir = root
        self.count = 0
        self.items_per_file = items_per_file
        self.files_per_dir = files_per_dir
        self.func = func
        self.indent = indent
        self.items = []

    def add(self, item):
        self.items.append(item)

        if len(self.items) == self.items_per_file:
            self.flush()

    def flush(self) -> Optional[str]:
        if not self.items:
            return None
        elif self.count + 1 == self.files_per_dir:
            # Too many files in directory: create a subdirectory
            self.dir = mkdtemp(dir=self.dir)
            os.chmod(self.dir, 0o775)
            self.count = 0

        fd, path = mkstemp(dir=self.dir)
        os.close(fd)
        os.chmod(path, 0o775)

        with open(path, "wt") as fh:
            json.dump(self.func(self.items), fh, indent=self.indent)

        self.count += 1
        self.items = []
        os.rename(path, path + ".json")
        return path + ".json"
