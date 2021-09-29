import bisect
import gzip
import multiprocessing as mp
import os
import pickle
import struct
import zlib
from typing import Callable, Optional, Sequence, Tuple

from .tempdir import TemporaryDirectory


class SimpleStore:
    def __init__(self, **kwargs):
        self._file = kwargs.get("file")
        self._tempdir = TemporaryDirectory(root=kwargs.get("tempdir"))
        self._fh = None

        if not self._file:
            self._file = self._tempdir.mktemp()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close(remove=True)

    def __del__(self):
        self.close(remove=True)

    def __iter__(self):
        self.close()

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
        return self._tempdir.size

    def add(self, item):
        if self._fh is None:
            self._fh = gzip.open(self._file, "wb", compresslevel=6)

        pickle.dump(item, self._fh)

    def close(self, remove: bool = False):
        if remove:
            self._tempdir.remove()

        if self._fh is None:
            return

        self._fh.close()
        self._fh = None


class Store:
    def __init__(self, file: str, mode: str = "r", **kwargs):
        self._file = file
        self._keys = kwargs.get("keys", [])
        self._tempdir = TemporaryDirectory(root=kwargs.get("tempdir"))
        self._bufmaxsize = kwargs.get("buffersize", 1000000)
        self._bufcursize = 0

        self._cache = {}
        self._files = []
        self._offsets = []

        if mode == "r":
            pass
        elif mode == "w":
            if not self._keys:
                raise ValueError(f"'keys' argument mandatory in write mode")

            os.makedirs(os.path.dirname(self._file), exist_ok=True)
            open(self._file, "w").close()

            for _ in self._keys:
                self._files.append(self._tempdir.mktemp())
        else:
            raise ValueError(f"invalid mode: '{mode}'")

    @property
    def file_keys(self):
        return self._keys

    @property
    def size(self) -> int:
        return self._tempdir.size

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def __getitem__(self, item):
        pass

    def __iter__(self):
        return self.keys()

    def keys(self):
        pass

    def values(self):
        pass

    def items(self):
        pass

    def range(self, start, stop=None):
        pass

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
        with open(self._file, "wb") as fh:
            # Header (empty for now)
            offset = fh.write(struct.pack("<Q", 0))

            # Body
            if workers > 1:
                ctx = mp.get_context(method="spawn")
                with ctx.Pool(workers - 1) as pool:
                    iterable = [(file, apply) for file in self._files]

                    for bytes_obj in pool.imap(self.load_items, iterable):
                        self._offsets.append(offset)
                        offset += fh.write(struct.pack("<L", len(bytes_obj)))
                        offset += fh.write(bytes_obj)
            else:
                for file in self._files:
                    bytes_obj = self.load_items((file, apply))
                    self._offsets.append(offset)
                    offset += fh.write(struct.pack("<L", len(bytes_obj)))
                    offset += fh.write(bytes_obj)

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

    @staticmethod
    def dump_items(fh, items):
        bytes_obj = zlib.compress(pickle.dumps(items))
        fh.write(struct.pack("<L", len(bytes_obj)))
        fh.write(bytes_obj)

    @staticmethod
    def load_items(args: Tuple[str, Optional[Callable]]) -> bytes:
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

    def close(self):
        pass