import argparse
import csv
import os
import shutil
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('analysis', nargs='?')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--working-directory')
    parser.add_argument('--rscript', default=shutil.which('Rscript'))
    args, extra = parser.parse_known_args()
    root = Path(__file__).resolve().parents[1]
    with (root / 'provenance/run_catalog.tsv').open(encoding='utf-8') as handle:
        catalog = list(csv.DictReader(handle, delimiter='\t'))
    if args.list:
        for item in catalog:
            print(f'{item["id"]}\t{item["analysis"]}\t{item["script"]}')
        return
    choices = {item['id']: item for item in catalog}
    if args.analysis not in choices:
        parser.error('Choose an analysis ID shown by --list.')
    item = choices[args.analysis]
    script = root / item['script']
    working = Path(args.working_directory).resolve() if args.working_directory else root / item['working_directory']
    print(f'Script: {script}')
    print(f'Working directory: {working}')
    if args.check:
        print(f'Script exists: {script.is_file()}')
        print(f'Working directory exists: {working.is_dir()}')
        return
    if script.suffix.lower() != '.r':
        parser.error('Run this Python preparation/layout script directly from the repository root.')
    if not args.rscript:
        parser.error('Rscript was not found. Add it to PATH or pass --rscript.')
    if not working.is_dir():
        parser.error('Stage the study inputs in the working directory before running this analysis.')
    env = os.environ.copy()
    env['DJK_PROJECT_ROOT'] = str(root)
    env['DJK16S_PROJECT_ROOT'] = str(root)
    env['DJK_SCRIPT_PATH'] = str(script)
    if item['working_directory'] != '.' or args.working_directory:
        env['DJK16S_ACTIVE_ROOT'] = str(working)
    else:
        env.pop('DJK16S_ACTIVE_ROOT', None)
    if extra and extra[0] == '--':
        extra = extra[1:]
    expression = 'source(Sys.getenv("DJK_SCRIPT_PATH"), chdir = FALSE, encoding = "UTF-8")'
    result = subprocess.run([args.rscript, '--vanilla', '-e', expression, *extra], cwd=working, env=env)
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
