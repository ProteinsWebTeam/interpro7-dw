import os
import shutil

LOADING_FILE = "loading"


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


from . import index, relationship, search
