import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('analysis', nargs='?')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--step')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--data-root', default=os.environ.get('DJK_PROJECT_ROOT'))
    parser.add_argument('--working-directory')
    parser.add_argument('--rscript', default=shutil.which('Rscript'))
    parser.add_argument('--bash', default=shutil.which('bash'))
    arguments = sys.argv[1:]
    if '--' in arguments:
        separator = arguments.index('--')
        extra = arguments[separator + 1:]
        arguments = arguments[:separator]
    else:
        extra = []
    args = parser.parse_args(arguments)
    root = Path(__file__).resolve().parents[2]
    with Path(__file__).with_name('workflow.tsv').open(encoding='utf-8') as handle:
        catalog = list(csv.DictReader(handle, delimiter='\t'))
    if args.list:
        for item in catalog:
            if args.step is None or item['id'].split('.')[0] == args.step.zfill(2):
                print(f'{item["id"]}\t{item["script"]}\t{item["role"]}')
        return
    choices = {item['id']: item for item in catalog}
    if args.analysis not in choices:
        parser.error('Choose an analysis ID shown by --list.')
    if not args.data_root:
        parser.error('Pass --data-root or set DJK_PROJECT_ROOT to the study data directory.')
    item = choices[args.analysis]
    data_root = Path(args.data_root).expanduser().resolve()
    script = root / item['script']
    working = Path(args.working_directory).expanduser().resolve() if args.working_directory else data_root / item['working_directory']
    print(f'Script: {script}', flush=True)
    print(f'Data root: {data_root}', flush=True)
    print(f'Working directory: {working}', flush=True)
    if args.check:
        print(f'Script exists: {script.is_file()}')
        print(f'Working directory exists: {working.is_dir()}')
        raise SystemExit(0 if script.is_file() and working.is_dir() else 1)
    if not script.is_file() or not working.is_dir():
        parser.error('The script and staged data working directory must exist.')
    env = os.environ.copy()
    env['DJK_CODE_ROOT'] = str(root)
    env['DJK_PROJECT_ROOT'] = str(data_root)
    env['DJK16S_PROJECT_ROOT'] = str(data_root)
    env['DJK_SCRIPT_PATH'] = str(script)
    if item['working_directory'] != '.' or args.working_directory:
        env['DJK16S_ACTIVE_ROOT'] = str(working)
    else:
        env.pop('DJK16S_ACTIVE_ROOT', None)
    suffix = script.suffix.lower()
    if suffix == '.r':
        if not args.rscript:
            parser.error('Rscript was not found. Add it to PATH or pass --rscript.')
        expression = '.libPaths(unique(c(file.path(Sys.getenv("DJK_PROJECT_ROOT"), "Rlib"), file.path(Sys.getenv("DJK_PROJECT_ROOT"), "99_运行环境", "R_library"), .libPaths()))); source(Sys.getenv("DJK_SCRIPT_PATH"), chdir = FALSE, encoding = "UTF-8")'
        command = [args.rscript, '--vanilla', '-e', expression, *extra]
    elif suffix == '.py':
        command = [sys.executable, str(script), *extra]
    elif suffix == '.sh':
        if not args.bash:
            parser.error('Bash was not found. Use Linux Bash or pass --bash.')
        command = [args.bash, str(script), *extra]
    else:
        parser.error('Unsupported script type.')
    result = subprocess.run(command, cwd=working, env=env)
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
