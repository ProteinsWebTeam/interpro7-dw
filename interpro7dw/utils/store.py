import bisect
import copy
import gzip
import multiprocessing as mp
import os
import pickle
import struct
import zlib
from tempfile import mkstemp
from typing import Any, Callable, Optional, Sequence, Tuple

from .tempdir import TemporaryDirectory


OptStr = Optional[str]


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
    with gzip.open(file, "wb") as fh:
        pickle.dump(data, fh)


def loadobj(file: str) -> Any:
    with gzip.open(file, "rb") as fh:
        return pickle.load(fh)


class SimpleStore:
    def __init__(self, file: OptStr = None, tempdir: OptStr = None):
        self._file = file
        self._is_tmp = False
        self._fh = None

        if not self._file:
            if tempdir:
                os.makedirs(tempdir, exist_ok=True)

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

            os.makedirs(os.path.dirname(self._file), exist_ok=True)
            open(self._file, "w").close()

            self._tempdir = TemporaryDirectory(root=kwargs.get("tempdir"))

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

    def merge(self, workers: int = 1, apply: Optional[Callable] = None):
        self.dump()

        with open(self._file, "wb") as fh:
            # Header (empty for now)
            offset = fh.write(struct.pack("<Q", 0))

            # Body
            for bytes_obj in self._merge(workers, apply):
                self._offsets.append(offset)
                offset += fh.write(struct.pack("<L", len(bytes_obj)))
                offset += fh.write(bytes_obj)

            # Footer
            pickle.dump((self._keys, self._offsets), fh)

            # Write footer offset in header
            fh.seek(0)
            fh.write(struct.pack("<Q", offset))

    def _merge(self, workers: int = 1, apply: Optional[Callable] = None):
        if workers > 1:
            ctx = mp.get_context(method="spawn")
            with ctx.Pool(workers - 1) as pool:
                iterable = [(file, apply) for file in self._files]
                yield from pool.imap(self.merge_items, iterable, chunksize=1000)
        else:
            for file in self._files:
                yield self.merge_items((file, apply))

    def dump(self):
        # Default to first file
        file = self._files[0]
        fh = open(file, "ab")

        keys = sorted(self._cache.keys())
        items = []
        for key in keys:
            # Find in which file store the values for the current key
            i = bisect.bisect_right(self._keys, key) - 1
            if i < 0:
                raise KeyError(key)

            if self._files[i] != file:
                # Different file: dump items if any
                if items:
                    self.dump_items(fh, items)
                    items.clear()

                # Open correct file
                fh.close()
                file = self._files[i]
                fh = open(file, "ab")

            items.append((key, self._cache.pop(key)))

        if items:
            self.dump_items(fh, items)

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
    def dump_items(fh, items):
        bytes_obj = zlib.compress(pickle.dumps(items))
        fh.write(struct.pack("<L", len(bytes_obj)))
        fh.write(bytes_obj)

    @staticmethod
    def merge_items(args: Tuple[str, Optional[Callable]]) -> bytes:
        file, apply = args

        data = {}
        with open(file, "rb") as fh:
            while True:
                bytes_obj = fh.read(4)
                if bytes_obj:
                    n_bytes, = struct.unpack("<L", bytes_obj)
                    items = pickle.loads(zlib.decompress(fh.read(n_bytes)))
                    for key, values in items:
                        try:
                            data[key] += values
                        except KeyError:
                            data[key] = values
                else:
                    break

        if apply is not None:
            for key, values in data.items():
                data[key] = apply(values)

        return zlib.compress(pickle.dumps(data))

    @staticmethod
    def get_first(values: Sequence):
        return values[0]

    @staticmethod
    def merge_dicts(values: Sequence[dict]) -> dict:
        dst = {}
        for value in values:
            copy_dict(value, dst)

        return dst
