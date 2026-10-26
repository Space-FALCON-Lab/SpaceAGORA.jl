#!/usr/bin/env python3
"""Prepare an offline GRAM loading environment from verified source and exact pins.

Requires Python 3.11+ and local Git objects. Does not invoke Julia or Pkg, fetch
objects, install dependencies, copy licensed assets, or alter either checkout.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tomllib

GRAM_UUID = 'b50455af-6a46-4eae-bf92-8039261dd674'
ROOT = Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def git(repo, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith('GIT_')}
    env.update(GIT_NO_LAZY_FETCH='1', GIT_OPTIONAL_LOCKS='0',
               GIT_TERMINAL_PROMPT='0', GIT_ALLOW_PROTOCOL='')
    result = subprocess.run(['git', '-C', str(repo), *args], env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise ValueError('Local Git read failed: ' + result.stderr.decode().strip())
    return result.stdout


def committed_file(root, revision, name):
    data = git(root, 'show', f'{revision}:{name}')
    if (root / name).read_bytes() != data:
        raise ValueError(f'{name} differs from the committed SpaceAGORA revision')
    return data


def unique_entry(entries, name):
    choices = entries.get(name, [])
    if len(choices) != 1:
        raise ValueError(f'Missing or ambiguous root manifest entry: {name}')
    return choices[0]


def dependency_closure(project, manifest):
    entries = manifest['deps']
    selected = {}
    queue = list(project['deps'].items())
    while queue:
        name, expected_uuid = queue.pop()
        item = unique_entry(entries, name)
        if item['uuid'] != expected_uuid:
            raise ValueError(f'Dependency UUID mismatch: {name}')
        if name in selected:
            continue
        if 'path' in item:
            raise ValueError(f'Path dependency needs separate preparation: {name}')
        selected[name] = item
        deps = item.get('deps', [])
        if isinstance(deps, dict):
            queue.extend(deps.items())
        else:
            queue.extend((n, unique_entry(entries, n)['uuid']) for n in deps)
        # Julia resolves list-form extension weakdeps through same-manifest
        # entries. Dictionary-form weakdeps already contain their UUIDs.
        weakdeps = item.get('weakdeps', {})
        if isinstance(weakdeps, list):
            queue.extend((n, unique_entry(entries, n)['uuid']) for n in weakdeps)
    return selected


def prepared_manifest(root_text, selected):
    manifest = tomllib.loads(root_text)
    if manifest.get('manifest_format') != '2.0':
        raise ValueError('A version 2.0 root manifest is required')
    blocks = re.split(r'(?=^\[\[deps\.)', root_text, flags=re.M)
    by_name = {}
    for block in blocks[1:]:
        match = re.match(r'\[\[deps\.([A-Za-z0-9_]+)\]\]', block)
        if not match or match[1] in by_name:
            raise ValueError('Unsupported or duplicate manifest stanza')
        by_name[match[1]] = block
    if not set(selected) <= set(by_name):
        raise ValueError('Selected dependency stanza is missing')
    # The root project_hash describes a different project. Retaining it would
    # claim a resolution of this wrapper project that has not been performed.
    header = ('# Offline loading metadata; exact SpaceAGORA dependency stanzas.\n'
              '# No resolve, instantiate, installation or native validation.\n'
              '# Root project_hash omitted: it does not describe this project.\n'
              f'julia_version = {json.dumps(manifest["julia_version"])}\n'
              'manifest_format = "2.0"\n\n')
    text = header + ''.join(by_name[n] for n in sorted(selected))
    if tomllib.loads(text)['deps'] != {n: [v] for n, v in selected.items()}:
        raise ValueError('Prepared manifest differs from the selected root entries')
    return text


def verify_source(source, objects, revision):
    files = {}
    tree = git(objects, 'ls-tree', '-rz', revision, 'Project.toml', 'src')
    for record in tree.split(b'\0'):
        if not record:
            continue
        meta, raw_name = record.split(b'\t', 1)
        mode, kind, oid = meta.decode().split()
        name = raw_name.decode()
        if kind != 'blob' or mode not in ('100644', '100755'):
            raise ValueError(f'Expected a regular source file: {name}')
        local = source / name
        if local.is_symlink() or not local.is_file():
            raise ValueError(f'Missing or linked source file: {name}')
        data = local.read_bytes()
        if data != git(objects, 'cat-file', 'blob', oid):
            raise ValueError(f'Source differs from configured GRAM pin: {name}')
        files[name] = {'git_blob': oid, 'sha256': digest(data)}
    actual = {'Project.toml'} | {
        f.relative_to(source).as_posix() for f in (source / 'src').rglob('*')
        if f.is_file() or f.is_symlink()}
    if set(files) != actual or 'src/GRAMSuite.jl' not in files:
        raise ValueError('Source inventory differs from the configured GRAM pin')
    return files


def prepare(root, source, objects, output, mode='regular'):
    root, source, objects = (Path(p).resolve() for p in (root, source, objects))
    output = Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise ValueError('Output must be a new directory')
    destination = output.resolve()
    if any(destination.is_relative_to(p) for p in (root, source)):
        raise ValueError('Output must be outside the SpaceAGORA and GRAM source checkouts')
    if not output.parent.is_dir():
        raise ValueError('Output parent must already exist')
    revision = git(root, 'rev-parse', 'HEAD').decode().strip()
    root_project = committed_file(root, revision, 'Project.toml')
    root_bytes = committed_file(root, revision, 'Manifest.toml')
    pin_bytes = committed_file(root, revision, '.github/gramsuite-revisions')
    pins = {}
    for line in pin_bytes.decode().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts = line.split()
        if len(parts) != 2 or parts[0] in pins or not re.fullmatch(r'[0-9a-f]{40}', parts[1]):
            raise ValueError('Invalid configured GRAM revisions')
        pins[parts[0]] = parts[1]
    if mode not in pins:
        raise ValueError(f'Missing configured GRAM mode: {mode}')
    pin = pins[mode]
    if mode == 'regular' and git(root, 'rev-parse', f'{revision}:data/GRAMSuite.jl').decode().strip() != pin:
        raise ValueError('Regular GRAM pin differs from the committed submodule gitlink')
    files = verify_source(source, objects, pin)
    project_bytes = (source / 'Project.toml').read_bytes()
    project = tomllib.loads(project_bytes.decode())
    if project.get('name') != 'GRAMSuite' or project.get('uuid') != GRAM_UUID:
        raise ValueError('Source does not identify the GRAMSuite package')
    selected = dependency_closure(project, tomllib.loads(root_bytes.decode()))
    text = prepared_manifest(root_bytes.decode(), selected)
    receipt = {
        'scope': 'Offline package-loading preparation only; native and worker readiness unverified',
        'space_revision': revision, 'gram_revision': pin, 'gram_mode': mode,
        'space_project': str(root), 'gram_source': str(source), 'gram_object_repository': str(objects),
        'source_files': files, 'root_project_sha256': digest(root_project),
        'root_manifest_sha256': digest(root_bytes), 'configured_pins_sha256': digest(pin_bytes),
        'prepared_project_sha256': digest(project_bytes), 'prepared_manifest_sha256': digest(text.encode()),
        'dependency_count': len(selected), 'dependencies': selected,
        'closure_policy': 'Strong deps and list-form weakdeps recursively; dictionary weakdeps unchanged',
        'source_link': str(source / 'src'), 'root_project_hash_omitted': True,
        'dependency_installation_performed': False,
    }
    # All source/metadata validation precedes the first output write.
    output.mkdir()
    (output / 'Project.toml').write_bytes(project_bytes)
    (output / 'Manifest.toml').write_text(text)
    (output / 'src').symlink_to(source / 'src', target_is_directory=True)
    (output / 'preparation.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--space-project', type=Path, default=ROOT)
    parser.add_argument('--gram-source', type=Path, required=True)
    parser.add_argument('--gram-object-repo', type=Path,
                        help='Local repository containing the pinned source objects; defaults to --gram-source')
    parser.add_argument('--mode', choices=('regular', 'dev'), default='regular')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = prepare(args.space_project, args.gram_source,
                         args.gram_object_repo or args.gram_source, args.output, args.mode)
    except (ValueError, OSError, KeyError, tomllib.TOMLDecodeError) as error:
        parser.exit(1, f'Preparation unavailable: {error}\n')
    print(f'Prepared {result["dependency_count"]} exact dependency entries for GRAM {result["gram_revision"]}')
    print(f'Receipt: {args.output / "preparation.json"}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
