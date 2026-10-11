"""Make a portable archive with a SHA256 inventory; omit runtime and logs."""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'exports' / 'Allen_Dataset_M1_V1.zip')
    args = parser.parse_args()
    excluded = {'.venv', '__pycache__', '.pytest_cache', '.cache', '.git', 'logs', 'exports'}
    large_raw_paths = [ROOT / 'input/merfish/cache', ROOT / 'input/rnaseq/raw']
    files = [p for p in ROOT.rglob('*') if p.is_file()
             and not excluded.intersection(p.relative_to(ROOT).parts)
             and not any(raw == p or raw in p.parents for raw in large_raw_paths)
             and not p.name.startswith('.env')
             and p.name != 'file_inventory.json'
             and p.resolve() != args.output.resolve()]
    inventory = [{'path': p.relative_to(ROOT).as_posix(),
                  'bytes': p.stat().st_size, 'sha256': digest(p)} for p in sorted(files)]
    invpath = ROOT / 'file_inventory.json'
    invpath.write_text(json.dumps(inventory, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(files) + [invpath]:
            archive.write(p, (Path('Allen Dataset_M1_V1') / p.relative_to(ROOT)).as_posix())
    with zipfile.ZipFile(args.output) as archive:
        assert archive.testzip() is None, 'Archive CRC verification failed'
        for item in inventory:
            h = hashlib.sha256()
            with archive.open('Allen Dataset_M1_V1/' + item['path']) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                    h.update(chunk)
            assert h.hexdigest() == item['sha256'], item['path']
    print(json.dumps({'archive': str(args.output), 'files': len(inventory),
                      'bytes': args.output.stat().st_size, 'sha256': digest(args.output)}, indent=2))


if __name__ == '__main__':
    main()
