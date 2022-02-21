import bisect
import copy
import os
import pickle
import shutil
import struct
from tempfile import mkdtemp
from typing import Callable, Optional


def copy_dict(src: dict, dst: dict, concat_or_incr: bool = False):
    for key, value in src.items():
        if key in dst:
            if isinstance(value, dict):
                copy_dict(value, dst[key], concat_or_incr)
            elif isinstance(value, (list, tuple)):
                dst[key] += value
            elif isinstance(value, set):
                dst[key] |= value
            elif isinstance(value, (int, float, str)) and concat_or_incr:
                dst[key] += value
            else:
                dst[key] = value
        else:
            dst[key] = copy.deepcopy(value)


class Directory:
    def __init__(self, root: Optional[str] = None,
                 tempdir: Optional[str] = None):
        if root:
            self.root = root
            self.keep = True
        else:
            if tempdir:
                os.makedirs(tempdir, mode=0o775, exist_ok=True)

            self.root = mkdtemp(dir=tempdir)
            self.keep = False

        os.chmod(self.root, 0o775)
        self.num_files = 0

    def get_size(self) -> int:
        if not self.root:
            return 0

        size = 0
        for root, dirs, files in os.walk(self.root):
            for name in files:
                size += os.path.getsize(os.path.join(root, name))

        return size

    def mktemp(self, suffix: str = "") -> str:
        self.num_files += 1
        filename = str(self.num_files).zfill(12) + suffix
        subdirs = [filename[:3], filename[3:6], filename[6:9]]
        dirpath = os.path.join(self.root, *subdirs)

        os.umask(0o002)
        os.makedirs(dirpath, mode=0o775, exist_ok=True)

        filepath = os.path.join(dirpath, filename)
        open(filepath, "w").close()
        os.chmod(filepath, 0o664)
        return filepath

    def remove(self):
        if self.keep or not self.root:
            return

        shutil.rmtree(self.root)
        self.root = None

    def __del__(self):
        self.remove()


class BasicStore:
    def __init__(self, file: str, mode: str = "r"):
        self.file = file
        self.fh = None

        if mode == "r":
            self.fh = open(self.file, "rb")
        elif mode in ("a", "w"):
            os.makedirs(os.path.dirname(self.file), mode=0o775, exist_ok=True)
            if mode == "w":
                self.fh = open(self.file, "wb")
        else:
            raise ValueError(f"invalid mode: '{mode}'")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        self.close()

        with open(self.file, "rb") as fh:
            while True:
                try:
                    obj = pickle.load(fh)
                except EOFError:
                    break
                else:
                    yield obj

    def write(self, obj):
        pickle.dump(obj, self.fh)

    def append(self, obj):
        with open(self.file, "ab") as fh:
            pickle.dump(obj, fh)

    def close(self):
        if self.fh:
            self.fh.close()


class KVStore:
    def __init__(self, file: str):
        self.file = file
        self._keys = []
        self.offsets = []
        self.offset = None
        self.cache = {}
        self.fh = open(self.file, "rb")

        offset, = struct.unpack("<Q", self.fh.read(8))
        self.fh.seek(offset)
        indices = pickle.load(self.fh)

        for key, offset in indices:
            self._keys.append(key)
            self.offsets.append(offset)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        return self.keys()

    def __getitem__(self, item):
        if item in self.cache:
            return self.cache[item]

        i = bisect.bisect_right(self._keys, item) - 1
        if i < 0:
            raise KeyError(item)

        try:
            offset = self.offsets[i]
        except IndexError:
            raise KeyError(item)

        self.load(offset)
        return self.cache[item]

    def get(self, item, default=None):
        try:
            return self[item]
        except KeyError:
            return default

    def range(self, start, stop=None):
        # Find in which bucket `start` is
        i = bisect.bisect_right(self._keys, start) - 1
        if i < 0:
            raise KeyError(start)

        # Load first bucket
        try:
            offset = self.offsets[i]
        except IndexError:
            raise KeyError(start)

        self.load(offset)

        # Yield items for first bucket
        for key in sorted(self.cache):
            if key < start:
                continue
            elif stop is None or key < stop:
                yield key, self.cache[key]
            else:
                return

        # If items in following buckets, load these buckets
        while True:
            i += 1

            try:
                offset = self.offsets[i]
            except IndexError:
                return

            self.load(offset)
            for key in sorted(self.cache):
                if key < start:
                    continue
                elif stop is None or key < stop:
                    yield key, self.cache[key]
                else:
                    return

    def close(self):
        self.fh.close()
        self.cache.clear()

    def get_keys(self):
        return self._keys

    def keys(self):
        self.offset = None
        for offset in self.offsets:
            self.load(offset)
            for key in sorted(self.cache):
                yield key

    def values(self):
        for key in self.keys():
            yield self.cache[key]

    def items(self):
        for key in self.keys():
            yield key, self.cache[key]

    def load(self, offset: int):
        if offset == self.offset:
            return

        self.fh.seek(offset)
        self.cache = pickle.load(self.fh)
        self.offset = offset


class KVStoreBuilder:
    def __init__(self, file: str, keys: list, tempdir: Optional[str] = None,
                 cachesize: int = 1000000):
        self.file = file
        self.keys = keys
        self.dir = Directory(tempdir=tempdir)
        self.max_cachesize = cachesize
        self.cachesize = 0
        self.cache = {}
        self.files = []
        self.indices = []

        for _ in self.keys:
            self.files.append(self.dir.mktemp())

        os.makedirs(os.path.dirname(self.file), mode=0o775, exist_ok=True)

        # Required if create store from sorted-data (without temp files)
        with open(self.file, "wb") as fh:
            # Header (empty for now)
            fh.write(struct.pack("<Q", 0))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.dir.remove()

    def __del__(self):
        self.dir.remove()

    def add(self, key, value):
        if key in self.cache:
            self.cache[key].append(value)
        else:
            self.cache[key] = [value]

        self.cachesize += 1
        if self.cachesize == self.max_cachesize:
            self.dump_to_tmp()

    def dump_to_tmp(self):
        file = self.files[0]
        fh = open(file, "ab")

        items = {}
        for key in sorted(self.cache.keys()):
            # Find in which file the values for the current key should be
            i = bisect.bisect_right(self.keys, key) - 1
            if i < 0:
                raise KeyError(key)

            if self.files[i] != file:
                # Different file: dump items and switch file
                if items:
                    pickle.dump(items, fh)
                    items.clear()

                fh.close()
                file = self.files[i]
                fh = open(file, "ab")

            items[key] = self.cache.pop(key)

        if items:
            pickle.dump(items, fh)

        fh.close()
        self.cachesize = 0

    def build(self, apply: Optional[Callable] = None):
        self.dump_to_tmp()

        with open(self.file, "wb") as fh:
            # Header (empty for now)
            fh.write(struct.pack("<Q", 0))

            # Body
            for key, file in zip(self.keys, self.files):
                self.indices.append((key, fh.tell()))

                data = self.load(file)
                if apply is not None:
                    for k, v in data.items():
                        data[k] = apply(v)

                pickle.dump(data, fh)

            # Footer
            offset = fh.tell()
            pickle.dump(self.indices, fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def append(self, key, value):
        self.cache[key] = value
        if len(self.cache) == self.max_cachesize:
            self.dump()

    def dump(self):
        if self.cache:
            with open(self.file, "ab") as fh:
                self.indices.append((min(self.cache.keys()), fh.tell()))
                pickle.dump(self.cache, fh)
                self.cache.clear()

    def get_size(self) -> int:
        return self.dir.get_size()

    def close(self):
        self.dump()

        with open(self.file, "ab") as fh:
            offset = fh.tell()
            pickle.dump(self.indices, fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    @staticmethod
    def load(file: str) -> dict:
        data = {}

        with open(file, "rb") as fh:
            while True:
                try:
                    obj = pickle.load(fh)
                except EOFError:
                    break
                else:
                    for key, values in obj.items():
                        if key in data:
                            data[key] += values
                        else:
                            data[key] = values

        return data

    @staticmethod
    def get_first(values: list):
        return values[0]
