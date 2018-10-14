#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import json
import logging
import os
import pickle
import shutil
import struct
import tempfile
import zlib
from typing import Iterable, Tuple


class Store(object):
    def __init__(self, filepath, verbose=False):
        self.filepath = filepath
        self.verbose = verbose
        self.mode = None
        self.keys = []
        self.offsets = []
        self.offset = 0
        self.data = {}

        if os.path.isfile(filepath):
            self.peek()

    def open(self, mode):
        if mode != self.mode:
            self.mode = mode

            if mode == 'w':
                with open(self.filepath, 'wb') as fh:
                    fh.write(struct.pack('<Q', 0))
                self.keys = []
                self.offsets = []
                self.offset = 8

    def add(self, data):
        if not data:
            return

        self.open('w')
        keys = sorted(data.keys())

        self.keys.append(keys[0].encode('utf-8'))
        self.offsets.append(self.offset)

        zstr = zlib.compress(json.dumps(data).encode('utf-8'))
        with open(self.filepath, 'ab') as fh:
            self.offset += fh.write(struct.pack('<I', len(zstr)) + zstr)

    def close(self):
        if self.mode == 'w':
            with open(self.filepath, 'rb+') as fh:
                fh.seek(0, 2)  # move at the end of the file
                fh.write(struct.pack('<I', len(self.keys)))

                for k, o in zip(self.keys, self.offsets):
                    fh.write(struct.pack('<IQ', len(k), o) + k)

                fh.seek(0)  # move at the start of the file
                fh.write(struct.pack('<Q', self.offset))

        self.mode = None
        self.keys = []
        self.offset = []
        self.offset = 0
        self.data = {}

    def peek(self):
        with open(self.filepath, 'rb') as fh:
            try:
                offset, = struct.unpack('<Q', fh.read(8))
            except struct.error:
                return False

            if not offset:
                return False

            fh.seek(offset)
            keys = []
            offsets = []

            try:
                n_keys, = struct.unpack('<I', fh.read(4))
                for _ in range(n_keys):
                    n_bytes, offset = struct.unpack('<IQ', fh.read(12))
                    keys.append(fh.read(n_bytes).decode('utf-8'))
                    offsets.append(offset)
            except struct.error:
                return False
            else:
                self.keys = keys
                self.offsets = offsets
                return True

    def get(self, k, default=None):
        if k in self.data:
            return self.data[k]
        elif not self.keys:
            raise RuntimeError("store at {} is empty".format(self.filepath))

        i = bisect.bisect_right(self.keys, k)
        if not i:
            return default

        offset = self.offsets[i-1]
        self.load(offset)
        return self.data.get(k, default)

    def load(self, offset, replace=True):
        if offset != self.offset:
            self.offset = offset

            if self.verbose:
                logging.info(
                    "{}: loading".format(os.path.basename(self.filepath))
                )

            with open(self.filepath, 'rb') as fh:
                fh.seek(offset)

                n_bytes, = struct.unpack('<I', fh.read(4))
                data = json.loads(zlib.decompress(fh.read(n_bytes)).decode('utf-8'))

                if replace:
                    self.data = data
                else:
                    self.data.update(data)

            if self.verbose:
                logging.info(
                    "{}: {} items loaded".format(
                        os.path.basename(self.filepath),
                        len(data)
                    )
                )

    def iter(self):
        if not self.offsets:
            raise RuntimeError("store at {} is empty".format(self.filepath))

        for offset in self.offsets:
            self.load(offset)

            for k in sorted(self.data):
                yield k, self.data[k]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()


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
            self.root = tempfile.mkdtemp(dir=kwargs.get("tmpdir"))
            self.n_aisles = kwargs.get("n_aisles")

            accessions = kwargs.get("accessions")
            accessions.sort()  # ensure the list is sorted

            data = []
            for i in range(0, len(accessions), self.n_aisles):
                fd, filepath = tempfile.mkstemp(dir=self.root)
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


class Bucket(object):
    def __init__(self, filepath, compress=False):
        self.filepath = filepath
        self.compress = compress
        self.keys = set()
        self.data = {}

    def add(self, key, _type, value):
        if key in self.data:
            if _type in self.data[key]:
                self.data[key][_type].add(value)
            else:
                self.data[key][_type] = {value}
        else:
            self.data[key] = {_type: value}
            self.keys.add(key)

    def dump(self):
        if self.data:
            if self.compress:
                s = zlib.compress(pickle.dumps(self.data))
            else:
                s = pickle.dumps(self.data)

            with open(self.filepath, "ab") as fh:
                fh.write(struct.pack("<I", len(s)) + s)

            self.data = {}

    def load():
        self.dump()
        data = {}

        with open(self.filepath, "rb") as fh:
            while True:
                try:
                    n_bytes, = struct.unpack("<I", fh.read(4))
                except struct.error:
                    break

                if self.compress:
                    chunk = pickle.loads(fh.read(n_bytes))
                else:
                    chunk = pickle.loads(zlib.decompress(fh.read(n_bytes)))

                for acc in chunk:
                    if acc in data:
                        for _type in chunk[acc]:
                            if _type in data[acc]:
                                data[acc][_type] |= chunk[acc][_type]
                            else:
                                data[acc][_type] = chunk[acc][_type]
                    else:
                        data[acc] = chunk[acc]

        return data

    def clean(self):
        os.remove(self.filepath)


class KVStore(object):
    def __init__(self, filepath, **kwargs):
        self.filepath = filepath
        self.bucket_size = kwargs.get("bucket_size", 1000)
        self.compress = kwargs.get("compress", False)
        self.tmpdir = kwargs.get("tmpdir")
        self.keys = {}
        self.buckets = []

    def add(self, key, _type, value):
        if key in self.keys:
            b = self.keys[key]
        else:
            try:
                b = self.buckets[-1]
            except IndexError:
                b = self.create_bucket()
            else:
                if len(b.keys) == self.bucket_size:
                    b = self.create_bucket()

            self.keys[key] = b

        b.add(key, _type, value)

    def create_bucket(self):
        fd, filepath = tempfile.mkstemp(dir=self.tmpdir)
        os.close(fd)

        b = Bucket(filepath, self.compress)
        self.buckets.append(b)
        return b

    def dump():
        for b in self.buckets:
            b.dump()

    def close():
        keys = {}
        offsets = {}
        with open(self.filepath, "wb") as fh:
            fh.write(struct.pack('<Q', 0))

            for b in self.buckets:
                data = b.load()
                b.clean()

                if self.compress:
                    s = zlib.compress(pickle.dumps(data))
                else:
                    s = pickle.dumps(data)

                o = fh.tell()
                fh.write(struct.pack("<Q", len(s)) + s)

                offsets[o] = []
                for k in data:
                    keys[k] = o
                    offsets[o].append(k)

            if self.compress:
                s = zlib.compress(pickle.dumps((keys, offsets)))
            else:
                s = pickle.dumps((keys, offsets))

            o = fh.tell()
            fh.write(struct.pack('<IQ', 1 if self.compress else 0, 0) + s)
            fh.seek(0)
            fh.write(struct.pack('<Q', o))
