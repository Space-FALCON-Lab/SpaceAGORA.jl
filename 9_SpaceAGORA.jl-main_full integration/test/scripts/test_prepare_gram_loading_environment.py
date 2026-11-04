#!/usr/bin/env python3
"""Offline preparation boundaries using tiny local Git/source fixtures.

Run with Python 3.11+; no Julia, package operations, licensed data or network.
"""
import copy
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location('prepare_gram', ROOT / 'scripts/prepare_gram_loading_environment.py')
helper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(helper)
A = '11111111-1111-1111-1111-111111111111'
B = '22222222-2222-2222-2222-222222222222'
W = '33333333-3333-3333-3333-333333333333'
D = '44444444-4444-4444-4444-444444444444'
ROOT_MANIFEST = f'''julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "not-the-wrapper-project"

[[deps.A]]
deps = ["B"]
uuid = "{A}"
version = "1.2.3"
weakdeps = ["W"]

    [deps.A.extensions]
    AWeakExt = "W"

[[deps.B]]
uuid = "{B}"
version = "4.5.6"

    [deps.B.weakdeps]
    Optional = "55555555-5555-5555-5555-555555555555"

[[deps.W]]
deps = ["D"]
uuid = "{W}"
version = "7.8.9"

[[deps.D]]
uuid = "{D}"
version = "10.11.12"
'''


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='gram-preparation-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name).resolve()
        self.root = self.base / 'space'
        self.source = self.base / 'gram'
        self.output = self.base / 'prepared'
        self.env = {k: v for k, v in os.environ.items() if not k.startswith('GIT_')}
        self.env.update(GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM='1',
                        GIT_TERMINAL_PROMPT='0', GIT_ALLOW_PROTOCOL='',
                        GIT_AUTHOR_NAME='Fixture', GIT_COMMITTER_NAME='Fixture',
                        GIT_AUTHOR_EMAIL='fixture@example.invalid', GIT_COMMITTER_EMAIL='fixture@example.invalid')
        for repo in (self.root, self.source):
            repo.mkdir()
            self.git(repo, 'init', '--quiet', '--template=')
        (self.source / 'src').mkdir()
        (self.source / 'src/GRAMSuite.jl').write_text('module GRAMSuite\nend\n')
        (self.source / 'Project.toml').write_text(
            f'name = "GRAMSuite"\nuuid = "{helper.GRAM_UUID}"\nversion = "0.1.0"\n[deps]\nA = "{A}"\n')
        self.commit(self.source)
        self.pin = self.git(self.source, 'rev-parse', 'HEAD').strip()
        (self.root / '.github').mkdir()
        (self.root / '.github/gramsuite-revisions').write_text(f'regular {self.pin}\ndev {self.pin}\n')
        (self.root / 'Project.toml').write_text('name = "SpaceAGORA"\n')
        (self.root / 'Manifest.toml').write_text(ROOT_MANIFEST)
        self.git(self.root, 'add', '.')
        self.git(self.root, 'update-index', '--add', '--cacheinfo', f'160000,{self.pin},data/GRAMSuite.jl')
        self.git(self.root, 'commit', '--quiet', '-m', 'fixture')

    def git(self, repo, *args):
        return subprocess.check_output(['git', '-C', str(repo), *args], env=self.env, stderr=subprocess.PIPE).decode()

    def commit(self, repo):
        self.git(repo, 'add', '.')
        self.git(repo, 'commit', '--quiet', '-m', 'fixture')

    def run_prepare(self, **kwargs):
        return helper.prepare(self.root, kwargs.get('source', self.source),
                              kwargs.get('objects', self.source),
                              kwargs.get('output', self.output), kwargs.get('mode', 'regular'))

    def rejected(self, reason, **kwargs):
        with self.assertRaisesRegex(ValueError, reason):
            self.run_prepare(**kwargs)
        self.assertFalse(self.output.exists())

    def test_exact_metadata_and_source_link_without_checkout_changes(self):
        before = {p: p.read_bytes() for repo in (self.root, self.source)
                  for p in repo.rglob('*') if p.is_file() and '.git' not in p.parts}
        receipt = self.run_prepare()
        prepared = tomllib.loads((self.output / 'Manifest.toml').read_text())
        self.assertEqual(set(prepared['deps']), {'A', 'B', 'W', 'D'})
        self.assertEqual(prepared['deps'], tomllib.loads(ROOT_MANIFEST)['deps'])
        self.assertNotIn('project_hash', prepared)
        self.assertEqual(receipt['gram_revision'], self.pin)
        self.assertEqual(receipt['dependency_count'], 4)
        self.assertEqual((self.output / 'src').resolve(), self.source / 'src')
        self.assertEqual((self.output / 'Project.toml').read_bytes(), (self.source / 'Project.toml').read_bytes())
        self.assertEqual(json.loads((self.output / 'preparation.json').read_text()), receipt)
        self.assertTrue(all(p.read_bytes() == data for p, data in before.items()))
        self.assertEqual(self.git(self.source, 'status', '--porcelain'), '')
        self.assertFalse(receipt['dependency_installation_performed'])

    def test_source_snapshot_requires_local_exact_pin_objects(self):
        snapshot = self.base / 'snapshot'
        snapshot.mkdir()
        (snapshot / 'src').mkdir()
        for name in ('Project.toml', 'src/GRAMSuite.jl'):
            (snapshot / name).write_bytes((self.source / name).read_bytes())
        result = self.run_prepare(source=snapshot)
        self.assertEqual(result['gram_revision'], self.pin)
        self.assertEqual((self.output / 'src').resolve(), snapshot / 'src')

    def test_modified_source_rejected_before_output(self):
        (self.source / 'src/GRAMSuite.jl').write_text('changed\n')
        self.rejected('Source differs')

    def test_extra_source_rejected(self):
        (self.source / 'src/injected.jl').write_text('injected\n')
        self.rejected('Source inventory differs')

    def test_missing_pin_objects_rejected_without_fetch(self):
        self.rejected('Local Git read failed', objects=self.root)

    def test_dirty_root_manifest_rejected(self):
        with (self.root / 'Manifest.toml').open('a') as stream:
            stream.write('# local metadata change\n')
        self.rejected('Manifest.toml differs')

    def test_dirty_root_project_rejected(self):
        with (self.root / 'Project.toml').open('a') as stream:
            stream.write('# local dependency change\n')
        self.rejected('Project.toml differs')

    def test_dirty_pin_configuration_rejected(self):
        with (self.root / '.github/gramsuite-revisions').open('a') as stream:
            stream.write('# local pin change\n')
        self.rejected('gramsuite-revisions differs')

    def test_gitlink_pin_disagreement_rejected(self):
        self.git(self.root, 'update-index', '--cacheinfo', f'160000,{"a" * 40},data/GRAMSuite.jl')
        self.git(self.root, 'commit', '--quiet', '-m', 'different link')
        self.rejected('pin differs.*gitlink')

    def test_existing_output_preserved(self):
        self.output.mkdir()
        sentinel = self.output / 'keep'
        sentinel.write_text('keep this')
        with self.assertRaisesRegex(ValueError, 'new directory'):
            self.run_prepare()
        self.assertEqual(sentinel.read_text(), 'keep this')

    def test_output_inside_checkout_rejected(self):
        self.rejected('outside', output=self.root / 'prepared')
        self.rejected('outside', output=self.source / 'prepared')

    def test_output_symlink_into_checkout_rejected(self):
        link = self.base / 'alias'
        link.symlink_to(self.root, target_is_directory=True)
        self.rejected('outside', output=link / 'prepared')

    def test_linked_source_file_rejected(self):
        saved = self.base / 'source.jl'
        saved.write_bytes((self.source / 'src/GRAMSuite.jl').read_bytes())
        (self.source / 'src/GRAMSuite.jl').unlink()
        (self.source / 'src/GRAMSuite.jl').symlink_to(saved)
        self.rejected('Missing or linked source')

    def test_missing_and_ambiguous_dependencies_rejected(self):
        for change, expected in [('missing', 'Missing or ambiguous'), ('duplicate', 'Missing or ambiguous'),
                                 ('uuid', 'UUID mismatch'), ('path', 'Path dependency')]:
            manifest = tomllib.loads(ROOT_MANIFEST)
            if change == 'missing':
                del manifest['deps']['A']
            elif change == 'duplicate':
                manifest['deps']['A'].append(copy.deepcopy(manifest['deps']['A'][0]))
            elif change == 'uuid':
                manifest['deps']['A'][0]['uuid'] = B
            else:
                manifest['deps']['A'][0]['path'] = '../elsewhere'
            with self.subTest(change=change), self.assertRaisesRegex(ValueError, expected):
                helper.dependency_closure({'deps': {'A': A}}, manifest)

    def test_list_weakdependency_missing_stanza_is_not_silently_dropped(self):
        manifest = tomllib.loads(ROOT_MANIFEST)
        del manifest['deps']['W']
        with self.assertRaisesRegex(ValueError, 'Missing or ambiguous.*W'):
            helper.dependency_closure({'deps': {'A': A}}, manifest)

    def test_cli_uses_documented_arguments(self):
        result = subprocess.run([os.sys.executable, str(ROOT / 'scripts/prepare_gram_loading_environment.py'),
                                 '--space-project', str(self.root), '--gram-source', str(self.source),
                                 '--output', str(self.output)], env=self.env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(self.pin, result.stdout)
        self.assertTrue((self.output / 'preparation.json').is_file())


if __name__ == '__main__':
    unittest.main()
