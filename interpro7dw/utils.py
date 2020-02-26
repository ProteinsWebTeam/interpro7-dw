# -*- coding: utf-8 -*-

import bisect
import copy
import multiprocessing as mp
import os
import pickle
import re
import shutil
import sqlite3
import struct
import tempfile
import zlib
from typing import Callable, Optional, Sequence, Tuple


class DirectoryTree(object):
    def __init__(self, root: Optional[str]=None, limit: int=1000):
        if root:
            os.makedirs(root, exist_ok=True)
            os.chmod(root, 0o775)

        self.root = tempfile.mkdtemp(dir=root)
        self.limit = limit
        self.cwd = self.root
        self.cnt = 0

    def mktemp(self, suffix=None, prefix=None) -> str:
        if self.cnt + 1 == self.limit:
            # Too many entries in the current directory: create subdirectory
            self.cwd = tempfile.mkdtemp(dir=self.cwd)
            self.cnt = 0
            os.chmod(self.cwd, 0o775)

        self.cnt += 1
        fd, path = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=self.cwd)
        os.close(fd)
        os.remove(path)
        os.chmod(path, 0o775)
        return path

    def remove(self):
        if not self.root:
            return

        shutil.rmtree(self.root)
        self.root = None


def deepupdate(input: dict, output: dict, replace: bool=True):
    for key, value in input.items():
        if key in output:
            if isinstance(value, dict):
                deepupdate(value, output[key], replace=replace)
            elif isinstance(value, (list, tuple)):
                output[key] += value
            elif isinstance(value, set):
                output[key] |= value
            elif replace:
                output[key] = value
            else:
                output[key] += value
        else:
            output[key] = copy.deepcopy(value)


class Bucket(object):
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.data = {}

        # Ensure the file always exists
        open(self.filepath, "wb").close()

        self.type = 0

    def __setitem__(self, key, value):
        self.data[key] = value
        self.type = 0

    def items(self):
        with open(self.filepath, "rb") as fh:
            while True:
                bytes_object = fh.read(4)
                if bytes_object:
                    n_bytes, = struct.unpack("<L", bytes_object)
                    data = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                    for key, value in data.items():
                        yield key, value
                else:
                    break

    @property
    def size(self) -> int:
        try:
            return os.path.getsize(self.filepath)
        except FileNotFoundError:
            return 0

    def append(self, key, value):
        try:
            self.data[key].append(value)
        except KeyError:
            self.data[key] = [value]
        finally:
            self.type = 1

    def add(self, key, value):
        try:
            self.data[key].add(value)
        except KeyError:
            self.data[key] = {value}
        finally:
            self.type = 2

    def update(self, key, value, replace: bool):
        if key in self.data:
            deepupdate(value, self.data[key], replace=replace)
        else:
            self.data[key] = copy.deepcopy(value)

        self.type = 3 if replace else 4

    def sync(self):
        if self.data:
            with open(self.filepath, "ab") as fh:
                bytes_object = zlib.compress(pickle.dumps(self.data))
                fh.write(struct.pack("<L", len(bytes_object)))
                fh.write(bytes_object)

            self.data = {}

    def merge(self) -> dict:
        self.sync()

        try:
            if self.type == 0:
                return self._merge()
            elif self.type == 1:
                return self._merge_list()
            elif self.type == 2:
                return self._merge_set()
            elif self.type == 3:
                return self._merge_dict(replace=True)
            elif self.type == 4:
                return self._merge_dict(replace=False)
            else:
                raise ValueError(self.type)
        finally:
            os.remove(self.filepath)

    def _merge(self) -> dict:
        return {key: value for key, value in self.items()}

    def _merge_list(self) -> dict:
        data = {}
        for key, value in self.items():
            try:
                data[key] += value
            except KeyError:
                data[key] = value

        return data

    def _merge_set(self) -> dict:
        data = {}
        for key, value in self.items():
            try:
                data[key] |= value
            except KeyError:
                data[key] = value

        return data

    def _merge_dict(self, replace: bool) -> dict:
        data = {}
        for key, value in self.items():
            if key in data:
                deepupdate(value, data[key], replace=replace)
            else:
                data[key] = value

        return data


class Store(object):
    def __init__(self, filepath: str, keys: Optional[Sequence]=None,
                 dir: Optional[str]=None):
        self.filepath = filepath

        if keys:
            # Writing mode
            self.dir = DirectoryTree(dir)
            self.fh = None
            self._keys = keys
            self.offsets = []
            self.buckets = [Bucket(self.dir.mktemp()) for _ in self._keys]
        else:
            # Reading mode
            self.dir = None
            self.fh = open(self.filepath, "rb")
            footer_offset, = struct.unpack("<Q", self.fh.read(8))
            self.fh.seek(footer_offset)
            self._keys, self.offsets = pickle.load(self.fh)
            self.buckets = []

        # Only used in reading mode
        self.offset = None
        self.data = {}

    @property
    def chunks(self):
        return self._keys

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, key):
        if key in self.data:
            return self.data[key]

        i = bisect.bisect_right(self._keys, key)
        if not i:
            raise KeyError(key)

        try:
            offset = self.offsets[i-1]
        except IndexError:
            raise KeyError(key)

        if self.offset == offset:
            """
            We already loaded data at this offset: 
            if the item is not in `self.data` then it's not in the store
            """
            raise KeyError(key)

        self._load(offset)
        self.offset = offset
        return self.data[key]

    def __iter__(self):
        return self.keys()

    def keys(self):
        for offset in self.offsets:
            self._load(offset)
            for key in sorted(self.data):
                yield key

    def values(self):
        for key in self.keys():
            yield self.data[key]

    def items(self):
        for key in self.keys():
            yield key, self.data[key]

    def __setitem__(self, key, value):
        self._get_bucket(key)[key] = value

    def _get_bucket(self, key) -> Bucket:
        i = bisect.bisect_right(self._keys, key)
        if i:
            return self.buckets[i-1]
        else:
            raise KeyError(key)

    def _load(self, offset: int):
        if self.fh is None:
            self.fh = open(self.filepath, "rb")

        self.fh.seek(offset)
        n_bytes, = struct.unpack("<L", self.fh.read(4))
        self.data = pickle.loads(zlib.decompress(self.fh.read(n_bytes)))

    def _merge_mp(self, fn: Optional[Callable], processes: int):
        offset = 0
        self.offsets = []

        with open(self.filepath, "wb") as fh, mp.Pool(processes-1) as pool:
            # Header (empty for now)
            offset += fh.write(struct.pack("<Q", 0))

            # Body
            iterable = [(bucket, fn) for bucket in self.buckets]
            for bytes_object in pool.imap(self._merge_bucket, iterable):
                self.offsets.append(offset)
                offset += fh.write(struct.pack("<L", len(bytes_object)))
                offset += fh.write(bytes_object)

            # Footer (index)
            pickle.dump((self._keys, self.offsets), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def _merge_sp(self, fn: Optional[Callable]):
        offset = 0
        self.offsets = []

        with open(self.filepath, "wb") as fh:
            # Header (empty for now)
            offset += fh.write(struct.pack("<Q", 0))

            # Body
            for bucket in self.buckets:
                self.offsets.append(offset)

                data = bucket.merge()

                if fn is not None:
                    self._dapply(data, fn)

                bytes_object = zlib.compress(pickle.dumps(data))
                offset += fh.write(struct.pack("<L", len(bytes_object)))
                offset += fh.write(bytes_object)

            # Footer (index)
            pickle.dump((self._keys, self.offsets), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def add(self, key, value):
        self._get_bucket(key).add(key, value)

    def append(self, key, value):
        self._get_bucket(key).append(key, value)

    def close(self):
        self.data.clear()

        if self.dir is not None:
            self.dir.remove()
            self.dir = None

        if self.fh is not None:
            self.fh.close()
            self.fh = None

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def merge(self, fn: Optional[Callable]=None, processes: int=1) -> int:
        self.sync()

        size = sum([b.size for b in self.buckets])
        if processes > 1:
            self._merge_mp(fn, processes)
        else:
            self._merge_sp(fn)

        return size

    def sync(self):
        for b in self.buckets:
            b.sync()

    def update(self, key, value, replace: bool):
        self._get_bucket(key).update(key, value, replace=replace)

    @staticmethod
    def dump_keys(keys: Sequence, output: str):
        with open(output, "wb") as fh:
            pickle.dump(keys, fh)

    @staticmethod
    def load_keys(filepath: str):
        with open(filepath, "rb") as fh:
            return pickle.load(fh)

    @staticmethod
    def _dapply(data: dict, fn: Callable):
        for key, value in data.items():
            data[key] = fn(value)

    @staticmethod
    def _merge_bucket(args: Tuple[Bucket, Optional[Callable]]) -> bytes:
        bucket, fn = args
        data = bucket.merge()

        if fn is not None:
            Store._dapply(data, fn)

        return zlib.compress(pickle.dumps(data))


class KVdb(object):
    def __init__(self, filepath: str, writeback: bool=False):
        self.filepath = filepath
        self.writeback = writeback
        self.con = sqlite3.connect(self.filepath)
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS data (
                id TEXT PRIMARY KEY NOT NULL,
                value TEXT NOT NULL
            )
            """
        )
        self.cache = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, key):
        try:
            return self.cache[key]
        except KeyError:
            pass

        sql = "SELECT value FROM data WHERE id = ?"
        row = self.con.execute(sql, (key,)).fetchone()

        if row is None:
            raise KeyError(key)

        value = pickle.loads(row[0])
        if self.writeback:
            self.cache[key] = value

        return value

    def __setitem__(self, key, value):
        if self.writeback:
            self.cache[key] = value
        else:
            sql = "INSERT OR REPLACE INTO data (id, value) VALUES (?, ?)"
            self.con.execute(sql, (key, pickle.dumps(value)))
            self.con.commit()

    def __iter__(self):
        return self.keys()

    def keys(self):
        for key, in self.con.execute("SELECT id FROM data ORDER BY id"):
            yield key

    def values(self):
        for obj, in self.con.execute("SELECT value FROM data ORDER BY id"):
            yield pickle.loads(obj)

    def items(self):
        for row in self.con.execute("SELECT id, value FROM data ORDER BY id"):
            yield row[0], pickle.loads(row[1])

    def close(self):
        if self.con is None:
            return

        self.sync()
        self.con.close()
        self.con = None

    def sync(self):
        if not self.cache:
            return

        sql = "INSERT OR REPLACE INTO data (id, value) VALUES (?, ?)"
        self.con.executemany(sql, ((key, pickle.dumps(value))
                                   for key, value in self.cache.items()))
        self.con.commit()
        self.cache = {}


def url2dict(url: str) -> dict:
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', url)

    if m is None:
        raise RuntimeError(f"invalid connection string: {url}")

    return dict(
        user=m.group(1),
        passwd=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        db=m.group(5)
    )


def datadump(filepath: str, data):
    with open(filepath, "wb") as fh:
        bytes_object = zlib.compress(pickle.dumps(data))
        fh.write(struct.pack("<L", len(bytes_object)))
        fh.write(bytes_object)


def dataload(filepath: str):
    with open(filepath, "rb") as fh:
        n_bytes, = struct.unpack("<L", fh.read(4))
        return pickle.loads(zlib.decompress(fh.read(n_bytes)))
