#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import json
import os
import shutil
import struct
import sys
import tempfile
import time
import zlib


class Store(object):
    def __init__(self, filepath):
        self.filepath = filepath
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

        i = bisect.bisect_right(self.keys, k)
        if not i:
            return default

        offset = self.offsets[i-1]
        self.load(offset)
        return self.data.get(k, default)

    def load(self, offset, replace=True):
        if offset != self.offset:
            self.offset = offset

            with open(self.filepath, 'rb') as fh:
                fh.seek(offset)

                n_bytes, = struct.unpack('<I', fh.read(4))
                data = json.loads(zlib.decompress(fh.read(n_bytes)).decode('utf-8'))

                if replace:
                    self.data = data
                else:
                    self.data.update(data)

    def iter(self):
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


class File(object):
    def __init__(self, path):
        self.path = path
        self.truncated = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def write(self, data):
        if not data:
            return 0
        elif self.truncated:
            mode = 'ab'
        else:
            mode = 'wb'
            self.truncated = True

        with open(self.path, mode) as fh:
            zstr = zlib.compress(json.dumps(data).encode('utf-8'))
            return fh.write(struct.pack('<I', len(zstr)) + zstr)

    def iter(self):
        with open(self.path, 'rb') as fh:
            while True:
                try:
                    n_bytes, = struct.unpack('<I', fh.read(4))
                except struct.error:
                    break
                else:
                    data = json.loads(zlib.decompress(fh.read(n_bytes)).decode('utf-8'))

                    for key in data:
                        yield key, data[key]


class XrefBucket(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.entries = {}
        self.n_xref = 0

    def add(self, accession, xref):
        if accession in self.entries:
            e = self.entries[accession]
        else:
            e = self.entries[accession] = {}

        # number of new cross-ref (all databases) for this entry
        n_new_xref = 0

        for dbname in xref:
            if dbname in e:
                n_xref = len(e[dbname])  # current number of cross-ref to `dbname` for this entry
                e[dbname] |= xref[dbname]
                n_new_xref += len(e[dbname]) - n_xref
            else:
                e[dbname] = xref[dbname]
                n_new_xref += len(xref[dbname])

        self.n_xref += n_new_xref
        return n_new_xref

    def dump(self):
        if self.n_xref:
            obj = {
                acc: {
                    dbname: list(self.entries[acc][dbname])
                    for dbname in self.entries[acc]
                }
                for acc in self.entries
            }

            n_dumped = self.n_xref
            self.entries = {}
            self.n_xref = 0

            with open(self.filepath, 'ab') as fh:
                zstr = zlib.compress(json.dumps(obj).encode('utf-8'))
                fh.write(struct.pack('<I', len(zstr)) + zstr)

            return n_dumped
        else:
            return 0

    def load(self):
        self.entries = {}
        with open(self.filepath, 'rb') as fh:
            while True:
                try:
                    n_bytes, = struct.unpack('<I', fh.read(4))
                except struct.error:
                    break
                else:
                    data = json.loads(zlib.decompress(fh.read(n_bytes)).decode('utf-8'))

                    for acc in data:
                        if acc in self.entries:
                            e = self.entries[acc]
                        else:
                            e = self.entries[acc] = {}

                        for dbname in data[acc]:
                            if dbname in e:
                                e[dbname] |= set(data[acc][dbname])
                            else:
                                e[dbname] = set(data[acc][dbname])

    def get(self, accession):
        if accession in self.entries:
            for dbname in self.entries[accession]:
                yield dbname, self.entries[accession][dbname]

    def free(self):
        self.entries = {}

    def __eq__(self, other):
        if isinstance(other, XrefBucket):
            return self.filepath == other.filepath
        else:
            return False


class Attic(object):
    def __init__(self, accessions, workdir=None, persist=False, max_xref=1000000):
        self.root = tempfile.mkdtemp(dir=workdir)
        self.accessions = accessions
        self.buckets = []
        self.persist = persist
        self.max_xref = max_xref
        self.bucket = None
        self.n_xref = 0

        for _ in range(len(self.accessions)):
            fd, filepath = tempfile.mkstemp(dir=self.root)
            os.close(fd)
            self.buckets.append(XrefBucket(filepath))

    def put(self, entries):
        for accession in entries:
            i = bisect.bisect_right(self.accessions, accession)
            if not i:
                return

            bucket = self.buckets[i - 1]
            self.n_xref += bucket.add(accession, entries[accession])

            if self.n_xref >= self.max_xref:
                # Xref limit reached: dump buckets

                '''
                Do not dump buckets that have less than half their "share" of xref.
                
                e.g.    if you `max_xref` is 1M and we have 10 buckets, each bucket can have 100k xref.
                        if one bucket has less than 50k xref, it's not worth dumping it.
                '''
                min_xref = self.max_xref / len(self.buckets) / 2

                for bucket in self.buckets:
                    if bucket.n_xref >= min_xref:
                        self.n_xref -= bucket.dump()

    def get(self, accession):
        i = bisect.bisect_right(self.accessions, accession)
        if not i:
            return []

        bucket = self.buckets[i - 1]
        if bucket != self.bucket:
            if self.bucket:
                self.bucket.free()

            self.bucket = bucket
            self.bucket.load()

        return self.bucket.get(accession)

    def close(self):
        for bucket in self.buckets:
            self.n_xref -= bucket.dump()

    def clean(self):
        shutil.rmtree(self.root)

    def __del__(self):
        if not self.persist and os.path.isdir(self.root):
            self.clean()


def test_read(proteins_f, *args, limit=2000000):
    proteins = Store(proteins_f)

    stores = [Store(arg) for arg in args]

    ts = time.time()
    cnt = 0
    for acc, protein in proteins.iter():
        for s in stores:
            _ = s.get(acc)

        cnt += 1
        if cnt == limit:
            break

    sys.stderr.write('{} ({:.0f} proteins/sec)\n'.format(cnt, cnt // (time.time() - ts)))

    proteins.close()
    for s in stores:
        s.close()
