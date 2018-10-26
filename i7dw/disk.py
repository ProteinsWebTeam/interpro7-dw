#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import json
import os
import pickle
import shutil
import struct
import zlib
from multiprocessing import Pool, Process, Queue
from tempfile import mkdtemp, mkstemp
from typing import Any, Callable, Generator, Iterable, Tuple, Union


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
    def __init__(self, filepath: str, keys: list=list(), tmpdir=None):
        self.filepath = filepath

        # To find chunk in filepath
        self.keys = keys
        self.offsets = []

        if tmpdir is not None:
            os.makedirs(tmpdir, exist_ok=True)

        if self.keys:
            self.dir = mkdtemp(dir=tmpdir)

            # Buckets when creating the file
            self.buckets = [self.create_bucket() for _ in self.keys]
        else:
            self.dir = None
            self.buckets = []

        # Type of values stored (None, default: overwrite any existing value)
        self.type = None

        # Variables used when reading the file
        self.fh = None
        self.items = {}
        self.offset = None

        self._dir = self.dir
        self.dir_limit = 1000
        self.dir_count = 0

        self.peek()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, key):
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

    def __setitem__(self, key, value):
        b = self.get_bucket(key)
        b[key] = value
        self.type = None

    def __iter__(self) -> Generator:
        for offset in self.offsets:
            self.load_chunk(offset)
            for key in sorted(self.items):
                yield key, self.items[key]

    def iter(self, maxsize: int=1) -> Generator:
        q = Queue(maxsize=maxsize)
        p = Process(target=self._iter, args=(self.filepath, self.offsets, q))
        p.start()

        while True:
            items = q.get()
            if items is None:
                break

            for key, value in items:
                yield key, value

    @staticmethod
    def _iter(filepath: str, offsets: list, queue: Queue):
        with open(filepath, "rb") as fh:
            for offset in offsets:
                fh.seek(offset)

                n_bytes, = struct.unpack("<L", fh.read(4))
                items = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                queue.put([(key, items[key]) for key in sorted(items)])

        queue.put(None)

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

        if self.dir:
            shutil.rmtree(self.dir)

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
                         [b.filepath for b in self.buckets]), fh)

        self.dir = None

    def reload(self):
        filepath = self.filepath + ".save"
        with open(filepath, "rb") as fh:
            self.dir, self.keys, files = pickle.load(fh)

        self.buckets = [Bucket(f) for f in files]
        os.remove(filepath)


class XrefAisle(object):
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.entries = {}
        self.n_xrefs = 0

    def add(self, accession: str, databases: dict) -> int:
        if accession in self.entries:
            e = self.entries[accession]
        else:
            e = self.entries[accession] = {}

        # Number of new cross-references
        n_new = 0

        for db in databases:
            if db in e:
                n = len(e[db])
                e[db] |= databases[db]  # expects a set
                n_new += len(e[db]) - n
            else:
                e[db] = databases[db]
                n_new += len(e[db])

        self.n_xrefs += n_new
        return n_new

    def dump(self) -> int:
        if not self.n_xrefs:
            return 0

        obj = {
            acc: {
                db: list(self.entries[acc][db])
                for db in self.entries[acc]
            }
            for acc in self.entries
        }

        n = self.n_xrefs
        self.entries = {}
        self.n_xrefs = 0

        with open(self.filepath, "ab") as fh:
            zstr = zlib.compress(json.dumps(obj).encode("utf-8"))
            fh.write(struct.pack("<I", len(zstr)) + zstr)

        return n

    def load(self):
        self.entries = {}
        with open(self.filepath, "rb") as fh:
            while True:
                try:
                    n_bytes, = struct.unpack("<I", fh.read(4))
                except struct.error:
                    break
                else:
                    data = json.loads(
                        zlib.decompress(fh.read(n_bytes)).decode("utf-8")
                    )

                    for acc in data:
                        if acc in self.entries:
                            e = self.entries[acc]
                        else:
                            e = self.entries[acc] = {}

                        for db in data[acc]:
                            if db in e:
                                e[db] |= set(data[acc][db])
                            else:
                                e[db] = set(data[acc][db])

    def get(self, accession: str) -> Iterable[Tuple[str, set]]:
        for db in self.entries.get(accession, {}):
            yield db, self.entries[accession][db]

    def free(self):
        self.entries = {}

    def __eq__(self, other):
        if isinstance(other, XrefAisle):
            return self.filepath == other.filepath
        else:
            return False


class XrefStore(object):
    def __init__(self, **kwargs):
        self.root = kwargs.get("root")
        self.max_xrefs = kwargs.get("max_xrefs", 1000000)
        self.accessions = []
        self.aisles = []
        self.aisle = None

        if self.root:
            # When already created (e.g. from child process)
            with open(os.path.join(self.root, "aisles.json"), "rt") as fh:
                for item in json.load(fh):
                    self.accessions.append(item["accession"])
                    self.aisles.append(XrefAisle(item["path"]))

            self.n_aisles = len(self.accessions)
        else:
            """
            When creating (from parent process)
            The following kwargs are required:
                - accessions
                - n_aisles
            """
            self.root = mkdtemp(dir=kwargs.get("tmpdir"))
            self.n_aisles = kwargs.get("n_aisles")

            accessions = kwargs.get("accessions")
            accessions.sort()  # ensure the list is sorted

            data = []
            for i in range(0, len(accessions), self.n_aisles):
                fd, filepath = mkstemp(dir=self.root)
                os.close(fd)
                self.accessions.append(accessions[i])
                self.aisles.append(XrefAisle(filepath))
                data.append({
                    "accession": accessions[i],
                    "path": filepath
                })

            with open(os.path.join(self.root, "aisles.json"), "wt") as fh:
                json.dump(data, fh)

        self.n_xrefs = 0

        """
        Do not dump aisles that have less
        than half their share of cross-refs

        e.g. if max_xrefs = 1M and n_aisles = 10,
            each aisle can have 100k cross-refs in memory
            so if one has less than 50k, it's not worth dumping it
        """
        self.min_xrefs = self.max_xrefs / self.n_aisles / 2

    def add(self, entries: dict):
        for accession in entries:
            i = bisect.bisect_right(self.accessions, accession)
            if not i:
                continue

            aisle = self.aisles[i-1]
            self.n_xrefs += aisle.add(accession, entries[accession])

            if self.n_xrefs >= self.max_xrefs:
                # Limit reached: dump cross-refs in memory to files
                for aisle in self.aisles:
                    if aisle.n_xrefs >= self.min_xrefs:
                        self.n_xrefs -= aisle.dump()

    def get(self, accession: str) -> Iterable[Tuple[str, set]]:
        i = bisect.bisect_right(self.accessions, accession)
        if not i:
            return []

        aisle = self.aisles[i-1]
        if aisle != self.aisle:
            if self.aisle:
                self.aisle.free()

            self.aisle = aisle
            self.aisle.load()

        return self.aisle.get(accession)

    def get_size(self) -> int:
        return sum([os.path.getsize(a.filepath) for a in self.aisles])

    def clean(self):
        for aisle in self.aisles:
            self.n_xrefs -= aisle.dump()

    def close(self):
        shutil.rmtree(self.root)
