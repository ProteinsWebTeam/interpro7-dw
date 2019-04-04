import json
import os
import shutil
import tempfile

LOADING_FILE = "loading"


class JsonFileOrganizer(object):
    def __init__(self, root: str, items_per_file: int=10000,
                 files_per_dir: int=1000):
        self.root = root  # must already exist
        self.dir = root
        self.count = 0
        self.items_per_file = items_per_file
        self.files_per_dir = files_per_dir
        self.items = []

    def add(self, item):
        self.items.append(item)

        if len(self.items) == self.items_per_file:
            self.flush()

    def flush(self):
        if not self.items:
            return
        elif self.count + 1 == self.files_per_dir:
            # Too many files in directory: create a subdirectory
            self.dir = tempfile.mkdtemp(dir=self.dir)
            self.count = 0

        fd, path = tempfile.mkstemp(dir=self.dir)
        os.close(fd)

        with open(path, "wt") as fh:
            json.dump(self.items, fh)

        self.count += 1
        self.items = []
        os.rename(path, path + ".json")


def init_dir(path: str):
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(path)
        open(os.path.join(path, LOADING_FILE), "w").close()


def is_ready(path: str):
    return not os.path.isfile(os.path.join(path, LOADING_FILE))


def set_ready(path: str):
    os.remove(os.path.join(path, LOADING_FILE))
