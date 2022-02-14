import re


def url2dict(uri: str) -> dict:
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', uri)

    if m is None:
        raise RuntimeError(f"invalid connection string: {uri}")

    return dict(
        user=m.group(1),
        passwd=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        db=m.group(5)
    )
