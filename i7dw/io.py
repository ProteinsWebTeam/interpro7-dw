# -*- coding: utf-8 -*-

import bisect
import copy
import os
import json
import pickle
import shutil
import sqlite3
import struct
import tempfile
import zlib
from multiprocessing import Pool, Queue
from typing import Any, Callable, Generator, Optional, Tuple, Union


def mktemp(prefix: Optional[str]=None, suffix: Optional[str]=None,
           dir: Optional[str]=None) -> str:
    if dir:
        os.makedirs(dir, exist_ok=True)

    fd, path = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=dir)
    os.close(fd)
    os.remove(path)

    return path


def serialize(value: dict) -> bytes:
    return pickle.dumps(value, pickle.HIGHEST_PROTOCOL)


def traverse(src: dict, dst: dict, replace: bool=True):
    for k, v in src.items():
        if k in dst:
            if isinstance(v, dict):
                traverse(v, dst[k], replace)
            elif isinstance(v, (list, tuple)):
                dst[k] += v
            elif isinstance(v, set):
                dst[k] |= v
            elif replace:
                dst[k] = v
            else:
                dst[k] += v
        else:
            dst[k] = copy.deepcopy(v)


class Bucket(object):
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.data = {}

        # Ensure the file always exists
        open(self.filepath, "wb").close()

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

    def update(self, key: Union[str, int], value: dict, replace: bool=True):
        if key in self.data:
            traverse(value, self.data[key], replace)
        else:
            self.data[key] = copy.deepcopy(value)

    def sync(self):
        if self.data:
            s = zlib.compress(serialize(self.data))
            self.data = {}
            with open(self.filepath, "ab") as fh:
                fh.write(struct.pack("<L", len(s)) + s)

    def merge_item(self) -> dict:
        return {k: v for k, v in self}

    def merge_list(self) -> dict:
        data = {}
        for k, v in self:
            if k in data:
                data[k] += v
            else:
                data[k] = v

        return data

    def merge_set(self) -> dict:
        data = {}
        for k, v in self:
            if k in data:
                data[k] |= v
            else:
                data[k] = v

        return data

    def merge_dict(self, replace: bool=True) -> dict:
        data = {}
        for k, v in self:
            if k in data:
                traverse(v, data[k], replace)
            else:
                data[k] = v

        return data

    def merge(self, merge_type: int) -> dict:
        self.sync()

        try:
            if merge_type == 0:
                return self.merge_item()
            elif merge_type == 1:
                return self.merge_list()
            elif merge_type == 2:
                return self.merge_set()
            elif merge_type in (3, 4):
                return self.merge_dict(replace=merge_type == 3)
            else:
                raise ValueError(merge_type)
        finally:
            os.remove(self.filepath)


class Store(object):
    def __init__(self, filepath: str, keys: Optional[list]=None,
                 tmpdir: Optional[str]=None, dir_limit: int=1000):
        self.filepath = filepath
        self.keys = keys

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
            3: update (replace numbers/strings)
            4: update (add numbers, concat strings)
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

            self._dir = self.dir = tempfile.mkdtemp(dir=tmpdir)

            # Create one bucket per key
            self.buckets = [Bucket(self._mktemp()) for _ in self.keys]
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
        return sum([bucket.size for bucket in self.buckets])

    def close(self):
        self.items = {}

        if self.dir:
            shutil.rmtree(self.dir)
            self.dir = None

        if self.fh is not None:
            self.fh.close()
            self.fh = None

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

    def append(self, key: Union[str, int], value: Any):
        self.get_bucket(key).append(key, value)
        self.type = 1

    def add(self, key: Union[str, int], value: Any):
        self.get_bucket(key).add(key, value)
        self.type = 2

    def update(self, key: Union[str, int], value: dict, replace: bool=True):
        self.get_bucket(key).update(key, value, replace)
        self.type = 3 if replace else 4

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

        return size

    def _merge_sp(self, func: Callable=None):
        pos = 0
        self.offsets = []

        with open(self.filepath, "wb") as fh:
            # Header (empty for now)
            pos += fh.write(struct.pack("<Q", 0))

            # Body
            for bucket in self.buckets:
                items = bucket.merge(self.type)

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

    def _merge_mp(self, processes: int, func: Callable=None):
        pos = 0
        self.offsets = []

        with open(self.filepath, "wb") as fh, Pool(processes) as pool:
            # Header (empty for now)
            pos += fh.write(struct.pack("<Q", 0))

            iterable = [(bucket, self.type, func) for bucket in self.buckets]
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
            self._dir = tempfile.mkdtemp(dir=self._dir)
            self.dir_count = 0

        self.dir_count += 1
        return mktemp(dir=self._dir)

    @staticmethod
    def _dapply(data: dict, func: Callable):
        for k, v in data.items():
            data[k] = func(v)

    @staticmethod
    def _merge_bucket(args: Tuple[Bucket, int, Optional[Callable]]) -> bytes:
        bucket, _type, func = args
        items = bucket.merge(_type)

        if func is not None:
            Store._dapply(items, func)

        return zlib.compress(serialize(items))

    @staticmethod
    def chunk_keys(keys, chunk_size: int) -> list:
        keys = sorted(keys)
        return [keys[i] for i in range(0, len(keys), chunk_size)]


class KVdb(object):
    def __init__(self, filepath: str, writeback: bool=False):
        self.filepath = filepath
        self.writeback = writeback
        self.con = sqlite3.connect(self.filepath)
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

    def __del__(self):
        self.close()

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

    def close(self):
        if self.con is None:
            return

        self.sync()
        self.con.close()
        self.con = None


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
            self.dir = tempfile.mkdtemp(dir=self.dir)
            os.chmod(self.dir, 0o775)
            self.count = 0

        path = mktemp(dir=self.dir)
        with open(path, "wt") as fh:
            json.dump(self.func(self.items), fh, indent=self.indent)

        os.chmod(path, 0o775)

        self.count += 1
        self.items = []
        os.rename(path, path + ".json")
        return path + ".json"
