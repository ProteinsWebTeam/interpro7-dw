# -*- coding: utf-8 -*-

import bisect
import copy
import heapq
import gzip
import multiprocessing as mp
import os
import pickle
import re
import shutil
import sqlite3
import struct
import zlib
from tempfile import mkdtemp, mkstemp
from typing import Callable, Iterable, Optional, Sequence, Tuple


class DirectoryTree:
    def __init__(self, root: Optional[str] = None, name: Optional[str] = None,
                 limit: int = 1000):
        if root:
            os.makedirs(root, exist_ok=True)

        if name:
            self.root = os.path.join(root, name)
            os.makedirs(self.root, exist_ok=True)
        else:
            self.root = mkdtemp(dir=root)

        os.chmod(self.root, 0o775)

        self.limit = limit
        self.cwd = self.root
        self.cnt = 0

    def mktemp(self, suffix: str = '', prefix: str = '') -> str:
        self.cnt += 1
        path = str(self.cnt).zfill(12)
        subdirs = [path[:3], path[3:6], path[6:9]]
        dirpath = os.path.join(self.root, *subdirs)

        os.umask(0o002)
        os.makedirs(dirpath, mode=0o775, exist_ok=True)

        filepath = os.path.join(dirpath, prefix + path + suffix)
        open(filepath, "w").close()
        os.chmod(filepath, 0o664)
        return filepath

    def remove(self):
        if not self.root:
            return

        shutil.rmtree(self.root)
        self.root = None

    @property
    def size(self) -> int:
        if not self.root:
            return 0

        size = 0
        for root, dirs, files in os.walk(self.root):
            for name in files:
                size += os.path.getsize(os.path.join(root, name))

        return size


def deepupdate(src: dict, dst: dict, replace: bool = True):
    for key, value in src.items():
        if key in dst:
            if isinstance(value, dict):
                deepupdate(value, dst[key], replace=replace)
            elif isinstance(value, (list, tuple)):
                dst[key] += value
            elif isinstance(value, set):
                dst[key] |= value
            elif replace:
                dst[key] = value
            else:
                dst[key] += value
        else:
            dst[key] = copy.deepcopy(value)


class Bucket:
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


class Store:
    def __init__(self, filepath: str, keys: Optional[Sequence] = None,
                 dir: Optional[str] = None):
        # Only used in reading mode
        self.offset = None
        self.data = {}
        self.num_items = 0

        if keys:
            # Writing mode
            self.dir = DirectoryTree(dir)
            self.filepath = filepath
            self.fh = None
            self._keys = keys
            self.offsets = []
            self.buckets = [Bucket(self.dir.mktemp()) for _ in self._keys]
        else:
            # Reading mode
            self.dir = None
            self.filepath = filepath
            self.fh = open(self.filepath, "rb")
            footer_offset, = struct.unpack("<Q", self.fh.read(8))
            self.fh.seek(footer_offset)
            self._keys, self.offsets, self.num_items = pickle.load(self.fh)
            self.buckets = []

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

    def __len__(self) -> int:
        return self.num_items

    def __setitem__(self, key, value):
        self._get_bucket(key)[key] = value

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

    def get_keys(self):
        return self._keys

    def keys(self):
        for offset in self.offsets:
            self._load(offset)
            for key in sorted(self.data):
                yield key

    def items(self):
        for key in self.keys():
            yield key, self.data[key]

    def merge(self, fn: Optional[Callable] = None, processes: int = 1) -> int:
        self.sync()

        size = sum([b.size for b in self.buckets])
        if processes > 1:
            self._merge_mp(fn, processes)
        else:
            self._merge_sp(fn)

        return size

    def range(self, start, stop: Optional):
        # Find in which bucket `start` is
        i = bisect.bisect_right(self._keys, start)
        if not i:
            raise KeyError(start)

        # Load first bucket
        offset = self.offsets[i-1]

        if offset != self.offset:
            self._load(offset)
            self.offset = offset

        # Yield items from the first bucket
        for key in sorted(self.data):
            if key < start:
                continue
            elif stop is None or key < stop:
                yield key, self.data[key]
            else:
                return

        # If items in following buckets, load these buckets
        while True:
            try:
                offset = self.offsets[i]
            except IndexError:
                return

            self._load(offset)
            self.offset = offset

            for key in sorted(self.data):
                if stop is None or key < stop:
                    yield key, self.data[key]
                else:
                    return

            i += 1

    def sync(self):
        for b in self.buckets:
            b.sync()

    def update(self, key, value, replace: bool):
        self._get_bucket(key).update(key, value, replace=replace)

    def values(self):
        for key in self.keys():
            yield self.data[key]

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
        self.num_items = 0

        with open(self.filepath, "wb") as fh:
            # Header (empty for now)
            offset += fh.write(struct.pack("<Q", 0))

            # Body
            ctx = mp.get_context(method="spawn")
            with ctx.Pool(processes-1) as pool:
                iterable = [(bucket, fn) for bucket in self.buckets]

                for bytes_obj, cnt in pool.imap(self._merge_bucket, iterable):
                    self.num_items += cnt
                    self.offsets.append(offset)
                    offset += fh.write(struct.pack("<L", len(bytes_obj)))
                    offset += fh.write(bytes_obj)

            # Footer
            pickle.dump((self._keys, self.offsets, self.num_items), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def _merge_sp(self, fn: Optional[Callable]):
        offset = 0
        self.offsets = []
        self.num_items = 0

        with open(self.filepath, "wb") as fh:
            # Header (empty for now)
            offset += fh.write(struct.pack("<Q", 0))

            # Body
            for bucket in self.buckets:
                data = bucket.merge()
                if fn is not None:
                    self._dapply(data, fn)

                self.num_items += len(data)
                self.offsets.append(offset)
                bytes_obj = zlib.compress(pickle.dumps(data))
                offset += fh.write(struct.pack("<L", len(bytes_obj)))
                offset += fh.write(bytes_obj)

            # Footer
            pickle.dump((self._keys, self.offsets, self.num_items), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    @staticmethod
    def chunk(keys: Iterable, chunk_size: int) -> Sequence:
        chunks = []
        for i, key in enumerate(sorted(keys)):
            if not i % chunk_size:
                chunks.append(key)

        return chunks

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
    def _merge_bucket(args: Tuple[Bucket, Optional[Callable]]) -> Tuple[bytes, int]:
        bucket, fn = args
        data = bucket.merge()

        if fn is not None:
            Store._dapply(data, fn)

        return zlib.compress(pickle.dumps(data)), len(data)


class KVdb:
    def __init__(self, filepath: str, writeback: bool = False):
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

    def __len__(self) -> int:
        return self.con.execute("SELECT COUNT(*) FROM data").fetchone()[0]

    def __delitem__(self, key):
        try:
            del self.cache[key]
        except KeyError:
            pass
        finally:
            self.con.execute("DELETE FROM data WHERE id = ?", (key,))
            self.con.commit()

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


def dumpobj(filepath: str, data):
    with open(filepath, "wb") as fh:
        bytes_object = zlib.compress(pickle.dumps(data))
        fh.write(struct.pack("<L", len(bytes_object)))
        fh.write(bytes_object)


def loadobj(filepath: str):
    with open(filepath, "rb") as fh:
        n_bytes, = struct.unpack("<L", fh.read(4))
        return pickle.loads(zlib.decompress(fh.read(n_bytes)))


class DumpFile:
    def __init__(self, path: str, compress: bool = False):
        self.path = path
        self.fh = None
        self.compresslevel = 6 if compress else 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        self.close()

        with gzip.open(self.path, "rb") as fh:
            while True:
                try:
                    obj = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield obj

    def close(self):
        if self.fh is None:
            return

        self.fh.close()
        self.fh = None

    def dump(self, obj):
        if self.fh is None:
            self.fh = gzip.open(self.path, "wb", self.compresslevel)

        pickle.dump(obj, self.fh)


def merge_dumps(files: Sequence[str], replace: bool = False):
    iterables = [DumpFile(path) for path in files]
    _key = None
    _values = None

    try:
        for key, values in heapq.merge(*iterables, key=lambda x: x[0]):
            if key == _key:
                deepupdate(values, _values, replace=replace)
            else:
                if _key is not None:
                    yield _key, _values

                _key = key
                _values = values

        if _key is not None:
            yield _key, _values
    finally:
        for df in iterables:
            df.close()
