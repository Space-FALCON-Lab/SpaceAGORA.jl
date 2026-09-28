#!/usr/bin/env python3
"""Test pinned GRAMSuite retrieval with independent local Git histories.

Run with Python 3 and Git; no third-party Python packages or network are used.
Pass --output-dir NEW_DIRECTORY to retain fixture repositories, logs and receipt.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

REPOSITORY = Path(__file__).resolve().parents[2]
HELPER = REPOSITORY / '.github/scripts/update_gramsuite_submodule.sh'
REVISIONS = REPOSITORY / '.github/gramsuite-revisions'


def test_retrieval(output):
    configured_pins = {}
    for line in REVISIONS.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parts = line.split()
        assert len(parts) == 2, 'Each revision row requires a mode and exact commit'
        mode, revision = parts
        assert mode in ('regular', 'dev') and mode not in configured_pins
        assert re.fullmatch(r'[0-9a-f]{40}', revision), 'Invalid configured pin'
        configured_pins[mode] = revision
    assert set(configured_pins) == {'regular', 'dev'}, 'Both pins are required'
    FIX = output / 'fixtures'
    FIX.mkdir()
    LOGS = output / 'test_logs'
    LOGS.mkdir()
    # Prevent inherited Git environment/configuration from reaching real repos,
    # credentials or remotes. URL rewrites route both production URLs to fixtures.
    ENV = {key: value for key, value in os.environ.items()
           if not key.startswith('GIT_') and key not in ('GH_TOKEN', 'GITHUB_TOKEN')}
    ENV.update(GIT_CONFIG_GLOBAL=str(FIX / 'gitconfig'), GIT_CONFIG_NOSYSTEM='1',
               GIT_TERMINAL_PROMPT='0', GIT_ALLOW_PROTOCOL='file', LC_ALL='C')
    RESULTS = []

    def run(*args, cwd=None, ok=True, env=None):
        p = subprocess.run(args, cwd=cwd, env=env or ENV, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if ok and p.returncode:
            raise AssertionError(f'{args}: {p.returncode}\n{p.stdout}{p.stderr}')
        return p

    def git(*args, cwd=None, ok=True):
        return run('git', *args, cwd=cwd, ok=ok)

    def commit(repo, message):
        git('add', '.', cwd=repo)
        git('-c', 'user.name=Pin fixture', '-c', 'user.email=pin@example.invalid',
            'commit', '--quiet', '-m', message, cwd=repo)
        return git('rev-parse', 'HEAD', cwd=repo).stdout.strip()

    pins = {}
    heads = {}
    for mode in ('regular', 'dev'):
        repo = FIX / (mode + '-remote')
        git('init', '--quiet', '--template=', '--initial-branch=main', str(repo))
        (repo / 'wrapper.txt').write_text(mode + ' accepted runtime\n')
        pins[mode] = commit(repo, mode + ' pinned')
        (repo / 'wrapper.txt').write_text(mode + ' later default HEAD\n')
        heads[mode] = commit(repo, mode + ' moved head')
        url = ('https://github.com/Space-FALCON-Lab/' +
               ('GRAMSuite.jl.git' if mode == 'regular' else 'dev-GRAMSuite.jl.git'))
        git('config', '--global', f'url.{repo.as_uri()}.insteadOf', url)
    git('config', '--global', 'protocol.file.allow', 'always')

    def checkout(name):
        repo = FIX / name
        git('init', '--quiet', '--template=', '--initial-branch=main', str(repo))
        (repo / '.github/scripts').mkdir(parents=True)
        shutil.copy2(HELPER, repo / '.github/scripts/update_gramsuite_submodule.sh')
        (repo / '.github/gramsuite-revisions').write_text(
            ''.join(f'{mode} {pins[mode]}\n' for mode in ('regular', 'dev')))
        (repo / '.gitmodules').write_text(
            '[submodule "GRAMSuite.jl"]\n'
            '\tpath = data/GRAMSuite.jl\n'
            '\turl = https://github.com/Space-FALCON-Lab/GRAMSuite.jl.git\n'
            '\tbranch = main\n')
        git('add', '.', cwd=repo)
        git('update-index', '--add', '--cacheinfo',
            f'160000,{pins["regular"]},data/GRAMSuite.jl', cwd=repo)
        git('-c', 'user.name=Pin fixture', '-c', 'user.email=pin@example.invalid',
            'commit', '--quiet', '-m', 'public gitlink with independent dev pin', cwd=repo)
        return repo

    def invoke(repo, mode, name, success=True, env=None):
        p = run('bash', str(repo / '.github/scripts/update_gramsuite_submodule.sh'),
                mode, cwd=FIX, ok=False, env=env)
        (LOGS / (name + '.log')).write_text(p.stdout + p.stderr)
        assert (p.returncode == 0) == success, (name, p.returncode, p.stdout, p.stderr)
        RESULTS.append({'name': name, 'pass': True, 'exit_code': p.returncode})
        return p

    def head(repo):
        return git('rev-parse', 'HEAD', cwd=repo / 'data/GRAMSuite.jl').stdout.strip()

    repo = checkout('fresh-regular')
    modules_before = (repo / '.gitmodules').read_bytes()
    invoke(repo, 'regular', 'fresh_regular_ignores_moved_default')
    assert head(repo) == pins['regular'] != heads['regular']
    assert git('symbolic-ref', '-q', 'HEAD', cwd=repo / 'data/GRAMSuite.jl', ok=False).returncode != 0
    assert (repo / '.gitmodules').read_bytes() == modules_before
    assert (repo / 'data/GRAMSuite.jl/wrapper.txt').read_text() == 'regular accepted runtime\n'

    repo_dev = checkout('fresh-dev')
    assert git('cat-file', '-e', pins['regular'], cwd=FIX / 'dev-remote', ok=False).returncode != 0
    invoke(repo_dev, 'dev', 'fresh_dev_without_public_gitlink_object')
    assert head(repo_dev) == pins['dev'] != heads['dev']
    assert (repo_dev / 'data/GRAMSuite.jl/wrapper.txt').read_text() == 'dev accepted runtime\n'
    assert git('rev-parse', '--git-dir', cwd=repo_dev / 'data/GRAMSuite.jl').stdout.strip().endswith('modules/GRAMSuite.jl')

    invoke(repo, 'dev', 'switch_regular_to_independent_dev')
    assert head(repo) == pins['dev']
    invoke(repo, 'regular', 'switch_dev_back_to_regular')
    assert head(repo) == pins['regular']
    invoke(repo, 'regular', 'repeat_same_pin')
    assert head(repo) == pins['regular']

    pin_file = repo / '.github/gramsuite-revisions'
    saved_pins = pin_file.read_text()
    pin_file.write_text(f'regular {pins["regular"]}\ndev {pins["regular"]}\n')
    failed = invoke(repo, 'dev', 'unavailable_exact_ref_no_fallback', success=False)
    assert head(repo) == pins['regular']
    assert 'not our ref' in failed.stderr or 'unadvertised object' in failed.stderr
    pin_file.write_text('regular PUBLIC_COMMIT_PENDING\ndev PRIVATE_COMMIT_PENDING\n')
    invoke(repo, 'regular', 'pending_pin_fails_before_fetch', success=False)
    assert head(repo) == pins['regular']
    pin_file.write_text(saved_pins + f'regular {pins["regular"]}\n')
    invoke(repo, 'regular', 'duplicate_pin_rejected', success=False)
    assert head(repo) == pins['regular']
    pin_file.write_text(saved_pins)
    invoke(repo, 'other', 'invalid_mode_rejected', success=False)

    (repo / 'data/GRAMSuite.jl/wrapper.txt').write_text('local edit must survive\n')
    invoke(repo, 'regular', 'dirty_same_pin_rejected', success=False)
    invoke(repo, 'dev', 'conflicting_local_edit_preserved', success=False)
    assert head(repo) == pins['regular']
    assert (repo / 'data/GRAMSuite.jl/wrapper.txt').read_text() == 'local edit must survive\n'
    git('checkout', '--', 'wrapper.txt', cwd=repo / 'data/GRAMSuite.jl')

    fake_token = 'fixture-token-this-is-not-a-credential'
    config_before = (FIX / 'gitconfig').read_bytes()
    p = invoke(repo, 'regular', 'token_not_persisted_or_logged',
               env=dict(ENV, GH_TOKEN=fake_token))
    assert fake_token not in p.stdout + p.stderr
    assert config_before == (FIX / 'gitconfig').read_bytes()
    for config in (repo / '.git/config', repo / '.git/modules/GRAMSuite.jl/config', repo / '.gitmodules'):
        assert fake_token not in config.read_text()
    assert not any('x-access-token:' in p.read_text() for p in LOGS.glob('*.log'))

    # Exercise the exact checked-in credential helper, without a network request.
    script_text = HELPER.read_text()
    helper_line = next(line for line in script_text.splitlines()
                       if "credential.https://github.com.helper=" in line)
    helper_arg = shlex.split(helper_line.rstrip().rstrip(chr(92)))[1]
    for host, expect_secret in [('github.com', True), ('example.invalid', False)]:
        result = subprocess.run(
            ['git', '-c', 'credential.helper=', '-c', helper_arg, 'credential', 'fill'],
            input=f'protocol=https\nhost={host}\n\n', text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env=dict(ENV, GH_TOKEN=fake_token))
        assert (fake_token in result.stdout) == expect_secret
        assert (result.returncode == 0) == expect_secret
        assert fake_token not in result.stderr
        RESULTS.append({'name': 'credential_scope_' + host,
                        'pass': True, 'exit_code': result.returncode})

    run('bash', '-n', str(HELPER))
    RESULTS.append({'name': 'bash_syntax', 'pass': True, 'exit_code': 0})
    receipt = {
        'completed_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'configured_wrapper_pins': configured_pins,
        'scope': 'Local Git fixture retrieval only; file transport is enforced. No Julia, native GRAM or real credentials.',
        'fixture_pins': pins, 'fixture_moved_default_heads': heads,
        'cases': RESULTS, 'all_passed': all(r['pass'] for r in RESULTS),
        'source_sha256': {
            str(p.relative_to(REPOSITORY)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [Path(__file__).resolve(), HELPER, REVISIONS]
        }
    }
    (output / 'TEST_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'cases': len(RESULTS), 'all_passed': receipt['all_passed'],
                      'configured_wrapper_pins': configured_pins}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        help='New directory in which to retain all test evidence')
    args = parser.parse_args()
    if args.output_dir:
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=False)
        test_retrieval(output)
    else:
        with tempfile.TemporaryDirectory(prefix='gramsuite-retrieval-') as work:
            test_retrieval(Path(work))


if __name__ == '__main__':
    main()
