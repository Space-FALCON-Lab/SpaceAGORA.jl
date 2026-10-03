"""Fail-closed export tests require only Python's standard library in PR CI."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'xval_fullarc_export.py'
spec = importlib.util.spec_from_file_location('xval_export', SCRIPT)
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)


class ExportProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.paths = []
        for target in ('gmat', 'stk'):
            run = self.root / f'{target}_committed'
            run.mkdir()
            self.paths.append(run)
            rows, cases = [], {}
            for name in sorted(export.SCENARIOS):
                case = run / name
                case.mkdir()
                (case / 'manifest.toml').write_text('fixture = true\n')
                (case / 'series.arrow').write_bytes(b'output identity fixture')
                ref = {'path': f'/retained/{target}/{name}.csv', 'sha256': 'a'*64}
                cases[name] = {'inputs': [ref], 'kernels': [ref], 'planetary_kernel': ref,
                               'solver_environment': {'mode': 'dp8'}, 'model': {'fixture': True},
                               'solver_retcode': 'Success',
                               'manifest_sha256': export.digest(case / 'manifest.toml'),
                               'series_sha256': export.digest(case / 'series.arrow')}
                rows.append({'scenario': name, 'target': target, 'variant': 'committed',
                             'reference_path': ref['path'], 'solver_retcode': 'Success',
                             'rms_m': 1.0, 'max_m': 2.0, 'n_points': 3, 'arc_end_s': 2.0,
                             'rms_first10k_m': 1.0, 'gravity_degree': 0, 'gravity_order': 0,
                             'gravity_file': 'fixture.csv', 'gm_override_m3s2': '', 'nbody_bodies': ''})
            with (run / 'results.csv').open('w', newline='') as stream:
                writer = csv.DictWriter(stream, rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            source = {'commit': 'b'*40, 'tree': 'c'*40, 'dirty': False}
            self.save(run, {'schema_version': 1, 'status': 'complete', 'target': target,
                            'variant': 'committed', 'source': source, 'source_end': source,
                            'controls': {}, 'scenarios': sorted(cases), 'cases': cases,
                            'reference_provenance': 'unverified_generation_settings',
                            'results_sha256': export.digest(run / 'results.csv')})

    def save(self, run, info):
        (run / 'run_info.json').write_text(json.dumps(info))

    def mutate(self, fn):
        run = self.paths[0]
        info = json.loads((run / 'run_info.json').read_text())
        fn(info)
        self.save(run, info)

    def rejected(self):
        with self.assertRaises((ValueError, KeyError, FileNotFoundError)):
            export.validate_campaign(self.root)

    def test_clean_exact_matrix(self):
        self.assertEqual(len(export.validate_campaign(self.root)), 2)

    def test_every_scientific_override_is_rejected(self):
        for key in ('XVAL_BASE', 'XVAL_GMAT_DIR', 'XVAL_J0_GM', 'XVAL_EARTH_FIELD_FILE', 'XVAL_MOON_FIELD_FILE',
                    'XVAL_MOON_FIELD_GM', 'XVAL_PLANETARY_KERNEL', 'XVAL_FRAME_TABLE_DIR',
                    'SPACEAGORA_SPICE_PCK_OVERRIDES', 'SPACEAGORA_GMAT_PARITY_SOLVER',
                    'SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH', 'XVAL_FUTURE_MODEL_OVERRIDE',
                    'SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT',
                    'SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT',
                    'SPACEAGORA_TELEMETRY_SOLVER_MAXITERS', 'SPACEAGORA_SOLVER_SAVE_EVERYSTEP',
                    'SPACEAGORA_NEW_UNKNOWN_CONTROL'):
            with self.subTest(key=key):
                self.mutate(lambda i: i.update(controls={key: 'diagnostic'}))
                self.rejected()

    def test_renamed_variant(self):
        self.mutate(lambda i: i.update(variant='as_basilisk'))
        self.rejected()

    def test_dirty_source_even_with_committed_label(self):
        self.mutate(lambda i: (i['source'].update(dirty=True), i['source_end'].update(dirty=True)))
        self.rejected()

    def test_mixed_source_revisions(self):
        self.mutate(lambda i: (i['source'].update(commit='d'*40), i['source_end'].update(commit='d'*40)))
        self.rejected()

    def test_failed_incomplete_missing_and_old_metadata(self):
        original = json.loads((self.paths[0] / 'run_info.json').read_text())
        for status in ('running', 'failed', ''):
            self.save(self.paths[0], dict(original, status=status))
            self.rejected()
        self.save(self.paths[0], dict(original, schema_version=0))
        self.rejected()
        (self.paths[0] / 'run_info.json').unlink()
        (self.paths[0] / 'run_info.toml').write_text('variant = "committed"')
        self.rejected()

    def test_stale_series_and_manifest(self):
        case = self.paths[0] / sorted(export.SCENARIOS)[0]
        for filename in ('series.arrow', 'manifest.toml'):
            path = case / filename
            before = path.read_bytes()
            path.write_bytes(b'changed')
            self.rejected()
            path.write_bytes(before)

    def test_changed_results_and_duplicate_coverage(self):
        path = self.paths[0] / 'results.csv'
        lines = path.read_text().splitlines(True)
        lines[-1] = lines[1]  # Still 24 rows; one scenario missing, another duplicated.
        path.write_text(''.join(lines))
        self.rejected()
        self.mutate(lambda i: i.update(results_sha256=export.digest(path)))
        self.rejected()

    def test_missing_input_and_unsuccessful_solver(self):
        name = sorted(export.SCENARIOS)[0]
        self.mutate(lambda i: i['cases'][name].update(inputs=[]))
        self.rejected()

    def test_loaded_planetary_kernel_identity_is_required(self):
        name = sorted(export.SCENARIOS)[0]
        self.mutate(lambda i: i["cases"][name].update(planetary_kernel={"path": "other.bsp", "sha256": "d"*64}))
        self.rejected()

    def test_failed_solver_even_if_record_and_rows_agree(self):
        path = self.paths[0] / 'results.csv'
        path.write_text(path.read_text().replace('Success', 'MaxIters'))
        def fail(info):
            info['results_sha256'] = export.digest(path)
            for case in info['cases'].values():
                case['solver_retcode'] = 'MaxIters'
        self.mutate(fail)
        self.rejected()

    def test_complete_export_preserves_configuration_and_provenance(self):
        try:
            import pandas
            import pyarrow as pa
            import pyarrow.feather as feather
        except ImportError:
            self.skipTest('Optional Arrow round trip; standard-library rejection tests still run')
        for run in self.paths:
            info = json.loads((run / 'run_info.json').read_text())
            for name, case in info['cases'].items():
                path = run / name / 'series.arrow'
                feather.write_feather(pa.table({'t_s': [0., 1., 2.], 'err_m': [0., 1., 2.]}), path)
                case['series_sha256'] = export.digest(path)
            self.save(run, info)
        paper = self.root / 'paper'
        (paper / 'data').mkdir(parents=True)
        (paper / 'figures' / 'validation_results').mkdir(parents=True)
        export.main(self.root, paper)
        with (paper / 'data' / 'cross_validation_rmse_fullarc_source.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 48)
        self.assertEqual(len({(r['body'], r['gravity'], r['third_body'], r['reference']) for r in rows}), 48)
        for row in rows:
            self.assertEqual(row['gravity_file'], 'fixture.csv')
            self.assertEqual(row['source_commit'], 'b'*40)
            self.assertEqual(row['reference_provenance'], 'unverified_generation_settings')
            self.assertEqual(len(row['run_record_sha256']), 64)

    def test_refusal_before_any_output(self):
        self.mutate(lambda i: i.update(status='running'))
        destination = self.root / 'paper'
        with self.assertRaises(ValueError):
            export.main(self.root, destination)
        self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()
