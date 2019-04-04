import json
import os
import tempfile


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
