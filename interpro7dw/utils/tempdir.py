import os
import shutil
from tempfile import mkdtemp
from typing import Optional


class TemporaryDirectory:
    def __init__(self, root: Optional[str] = None, keep: bool = False):
        self.keep = keep
        if root:
            os.makedirs(root, exist_ok=True)

        self.root = mkdtemp(dir=root)
        os.chmod(self.root, 0o775)

        self.num_files = 0

    @property
    def size(self) -> int:
        if not self.root:
            return 0

        size = 0
        for root, dirs, files in os.walk(self.root):
            for name in files:
                size += os.path.getsize(os.path.join(root, name))

        return size

    def mktemp(self) -> str:
        self.num_files += 1
        filename = str(self.num_files).zfill(12)
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
