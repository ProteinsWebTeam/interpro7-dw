import bisect
import copy
import gzip
import multiprocessing as mp
import os
import pickle
import shutil
import struct
import zlib
from tempfile import mkdtemp, mkstemp
from typing import Any, Callable, Optional, Sequence, Tuple


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


def dumpobj(data: Any, file: str):
    os.makedirs(os.path.dirname(file), mode=0o775, exist_ok=True)
    with gzip.open(file, "wb") as fh:
        pickle.dump(data, fh)


def loadobj(file: str) -> Any:
    with gzip.open(file, "rb") as fh:
        return pickle.load(fh)


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

    def size(self) -> int:
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


class SimpleStore:
    def __init__(self, file: Optional[str] = None,
                 tempdir: Optional[str] = None):
        self._file = file
        self._is_tmp = False
        self._fh = None

        if not self._file:
            if tempdir:
                os.makedirs(tempdir, mode=0o775, exist_ok=True)

            fd, self._file = mkstemp(dir=tempdir)
            os.close(fd)
            self._is_tmp = True

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __iter__(self):
        if self._fh:
            self._fh.close()
            self._fh = None

        with gzip.open(self._file, "rb") as fh:
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
            return os.path.getsize(self._file)
        except FileNotFoundError:
            return 0

    def add(self, item, keep_open: bool = True):
        if keep_open:
            if not self._fh:
                os.makedirs(os.path.dirname(self._file), mode=0o775,
                            exist_ok=True)
                self._fh = gzip.open(self._file, "wb", compresslevel=6)

            pickle.dump(item, self._fh)
        else:
            with gzip.open(self._file, "ab", compresslevel=6) as fh:
                pickle.dump(item, fh)

    def close(self):
        if self._fh:
            self._fh.close()
            self._fh = None

        if self._is_tmp and os.path.isfile(self._file):
            os.remove(self._file)


class Store:
    def __init__(self, file: str, mode: str = "r", **kwargs):
        self._file = file
        self._keys = kwargs.get("keys", [])
        self._bufmaxsize = kwargs.get("buffersize", 1000000)
        self._bufcursize = 0

        self._cache = {}
        self._files = []
        self._offsets = []
        self._tempdir = None
        self._fh = self._offset = None

        if mode == "r":
            self._keys, self._offsets = self.load_footer(self._file)
        elif mode == "w":
            if not self._keys:
                raise ValueError(f"'keys' argument mandatory in write mode")

            os.makedirs(os.path.dirname(self._file), mode=0o775, exist_ok=True)
            open(self._file, "w").close()

            self._tempdir = Directory(tempdir=kwargs.get("tempdir"))

            for _ in self._keys:
                self._files.append(self._tempdir.mktemp())
        else:
            raise ValueError(f"invalid mode: '{mode}'")

    @property
    def file_keys(self):
        return self._keys

    @property
    def size(self) -> int:
        return self._tempdir.size if self._tempdir else 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, item):
        if item in self._cache:
            return self._cache[item]

        i = bisect.bisect_right(self._keys, item) - 1
        if i < 0:
            raise KeyError(item)

        try:
            offset = self._offsets[i]
        except IndexError:
            raise KeyError(item)

        if self.load(offset):
            return self._cache[item]

        # Data at this offset already loaded
        raise KeyError(item)

    def __iter__(self):
        return self.keys()

    def close(self):
        if self._tempdir:
            self._tempdir.remove()
            self._tempdir = None

        if self._fh is None:
            return

        self._fh.close()
        self._fh = None

    def keys(self):
        self._offset = None
        for offset in self._offsets:
            self.load(offset)
            for key in sorted(self._cache):
                yield key

    def values(self):
        for key in self.keys():
            yield self._cache[key]

    def items(self):
        for key in self.keys():
            yield key, self._cache[key]

    def range(self, start, stop=None):
        # Find in which bucket `start` is
        i = bisect.bisect_right(self._keys, start) - 1
        if i < 0:
            raise KeyError(start)

        # Load first bucket
        try:
            offset = self._offsets[i]
        except IndexError:
            raise KeyError(start)

        self.load(offset)

        # Yield items for first bucket
        for key in sorted(self._cache):
            if key < start:
                continue
            elif stop is None or key < stop:
                yield key, self._cache[key]
            else:
                return

        # If items in following buckets, load these buckets
        while True:
            try:
                offset = self._offsets[i]
            except IndexError:
                return

            self.load(offset)
            for key in sorted(self._cache):
                if stop is None or key < stop:
                    yield key, self._cache[key]
                else:
                    return

            i += 1

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def add(self, key, value):
        try:
            self._cache[key].append(value)
        except KeyError:
            self._cache[key] = [value]

        self._bufcursize += 1
        if self._bufcursize == self._bufmaxsize:
            self.dump()

    def _merge_mp(self, fh, offset: int, workers: int,
                  apply: Optional[Callable]) -> int:
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=workers) as pool:
            iterable = [(file, apply) for file in self._files]

            for file in pool.imap(self.merge_items, iterable, chunksize=1000):
                self._offsets.append(offset)

                with open(file, "rb") as ifh:
                    offset += fh.write(ifh.read())

        return offset

    def _merge_sp(self, fh, offset: int, apply: Optional[Callable]) -> int:
        for file in self._files:
            self._offsets.append(offset)
            data = self.load_items(file)

            if apply is not None:
                for key, values in data.items():
                    data[key] = apply(values)

            bytes_obj = zlib.compress(pickle.dumps(data))
            offset += fh.write(struct.pack("<L", len(bytes_obj)))
            offset += fh.write(bytes_obj)

        return offset

    def merge(self, apply: Optional[Callable] = None, workers: int = 1):
        self.dump()

        with open(self._file, "wb") as fh:
            # Header (empty for now)
            offset = fh.write(struct.pack("<Q", 0))

            # Body
            if workers > 1:
                offset = self._merge_mp(fh, offset, workers - 1, apply)
            else:
                offset = self._merge_sp(fh, offset, apply)

            # Footer
            pickle.dump((self._keys, self._offsets), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def dump(self):
        # Default to first file
        file = self._files[0]
        fh = open(file, "ab")

        keys = sorted(self._cache.keys())
        items = {}
        for key in keys:
            # Find in which file store the values for the current key
            i = bisect.bisect_right(self._keys, key) - 1
            if i < 0:
                raise KeyError(key)

            if self._files[i] != file:
                # Different file: dump items if any
                if items:
                    bytes_obj = zlib.compress(pickle.dumps(items))
                    fh.write(struct.pack("<L", len(bytes_obj)))
                    fh.write(bytes_obj)
                    items.clear()

                # Open correct file
                fh.close()
                file = self._files[i]
                fh = open(file, "ab")

            items[key] = self._cache.pop(key)

        if items:
            bytes_obj = zlib.compress(pickle.dumps(items))
            fh.write(struct.pack("<L", len(bytes_obj)))
            fh.write(bytes_obj)

        fh.close()
        self._bufcursize = 0

    def load(self, offset: int) -> bool:
        if self._fh is None:
            self._fh = open(self._file, "rb")

        if offset == self._offset:
            return False  # Already loaded

        self._fh.seek(offset)
        n_bytes, = struct.unpack("<L", self._fh.read(4))
        self._cache = pickle.loads(zlib.decompress(self._fh.read(n_bytes)))
        self._offset = offset
        return True

    @staticmethod
    def load_footer(file: str):
        with open(file, "rb") as fh:
            footer_offset, = struct.unpack("<Q", fh.read(8))

            fh.seek(footer_offset)
            keys, offsets = pickle.load(fh)

        return keys, offsets

    @staticmethod
    def load_items(file: str) -> dict:
        data = {}
        with open(file, "rb") as fh:
            while True:
                bytes_obj = fh.read(4)
                if bytes_obj:
                    n_bytes, = struct.unpack("<L", bytes_obj)
                    chunk = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                    for key, values in chunk.items():
                        try:
                            data[key] += values
                        except KeyError:
                            data[key] = values
                else:
                    break

        return data

    @staticmethod
    def merge_items(args: Tuple[str, Optional[Callable]]):
        file, apply = args

        data = Store.load_items(file)

        if apply is not None:
            for key, values in data.items():
                data[key] = apply(values)

        with open(file, "wb") as fh:
            bytes_obj = zlib.compress(pickle.dumps(data))
            fh.write(struct.pack("<L", len(bytes_obj)))
            fh.write(bytes_obj)

        return file

    @staticmethod
    def get_first(values: Sequence):
        return values[0]


class StoreDirectLoader:
    def __init__(self, file: str, chunksize: int = 10000):
        self.fh = open(file, "wb")
        self.fh.write(struct.pack("<Q", 0))
        self.chunksize = chunksize
        self.cache = {}
        self.keys = []
        self.offsets = []

    def add(self, key, value):
        self.cache[key] = value
        if len(self.cache) == self.chunksize:
            self.dump()

    def dump(self):
        if not self.cache:
            return

        self.keys.append(min(self.cache.keys()))
        self.offsets.append(self.fh.tell())

        bytes_obj = zlib.compress(pickle.dumps(self.cache))
        self.fh.write(struct.pack("<L", len(bytes_obj)))
        self.fh.write(bytes_obj)

        self.cache.clear()

    def close(self):
        self.dump()

        footer_offset = self.fh.tell()
        pickle.dump((self.keys, self.offsets), self.fh)

        self.fh.seek(0)
        self.fh.write(struct.pack("<Q", footer_offset))
        self.fh.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.fh.close()

    def __del__(self):
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

        # If building directly from sorted data (no temp files)
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
        with open(self.file, "ab") as fh:
            self.indices.append((min(self.cache.keys()), fh.tell()))
            pickle.dump(self.cache, fh)
            self.cache.clear()

    def size(self) -> int:
        return self.dir.size()

    def close(self, direct: bool = False, apply: Optional[Callable] = None):
        if direct:
            self.dump()

            with open(self.file, "ab") as fh:
                offset = fh.tell()
                pickle.dump(self.indices, fh)

                # Write footer offset in header
                fh.seek(0)
                fh.write(struct.pack("<Q", offset))
        else:
            self.build(apply=apply)

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
