#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import os
import pickle
import shutil
import struct
import zlib
from multiprocessing import Process, Pool, Queue
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

    def items(self, processes: int=1) -> Callable:
        if processes > 1:
            return self.iter(processes-1)
        else:
            return self

    def iter(self, processes: int=1) -> Generator:
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
