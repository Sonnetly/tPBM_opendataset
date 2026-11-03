"""Extract Microsoft's Arial core fonts into a local cache; do not commit fonts."""
from pathlib import Path
import argparse
import hashlib
import json
import urllib.request
import libarchive

URL = 'https://downloads.sourceforge.net/corefonts/arial32.exe'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, default=Path('/tmp/ros_arial'))
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    data = urllib.request.urlopen(URL, timeout=120).read()
    exe = args.directory / 'arial32.exe'
    exe.write_bytes(data)
    hashes = {}
    with libarchive.file_reader(str(exe)) as archive:
        for entry in archive:
            name = Path(entry.pathname).name
            content = b''.join(entry.get_blocks())
            if name.lower().endswith('.ttf'):
                (args.directory / name).write_bytes(content)
                hashes[name] = hashlib.sha256(content).hexdigest()
    print(json.dumps({'source_url': URL, 'installer_sha256': hashlib.sha256(data).hexdigest(), 'fonts': hashes}, indent=2))

if __name__ == '__main__': main()
