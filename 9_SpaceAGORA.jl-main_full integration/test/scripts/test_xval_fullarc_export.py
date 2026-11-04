"""Fail-closed export tests need only Python 3.11's standard library in PR CI."""
import copy
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
                (case / 'series.arrow').write_bytes(b'output identity fixture')
                ref = {'path': f'/retained/{target}/{name}.csv', 'sha256': 'a'*64}
                body, grav, third_body = name.split('_')
                degree = {'j0': 0, 'j2': 2, 'j50': 50}[grav]
                bodies = {'earth': ['sun', 'moon'], 'mars': ['sun'],
                          'venus': ['sun'], 'moon': ['sun', 'earth']}[body] if third_body == 'tbtrue' else []
                field = f'data/Gravity_harmonics_data/{body}.csv'
                model = {'name': name, 'planet': body, 'gravity_model': 'inverse_squared',
                         'gravity_harmonics_degree': degree, 'gravity_harmonics_order': degree,
                         'gravity_harmonics_file': field, 'nbody_bodies': bodies,
                         'orbit_altitude_mode': 'oblate',
                         'initial_time': {'year': 2026, 'month': 1, 'day': 1,
                                          'hour': 12, 'minute': 0, 'second': 0.0},
                         'kind': 'time_aligned_state', 'comparison_mode': 'time_aligned_state',
                         'events': ['altitude_time', 'state_x_time', 'state_y_time', 'state_z_time'],
                         'telemetry': f'/retained/run/{name}/trajectory.csv',
                         'telemetry_columns': {'time': 'time_s', 'altitude': 'altitude_km',
                                               'x': 'x_km', 'y': 'y_km', 'z': 'z_km'},
                         'max_points_quick': 1_000_000_000, 'max_points_full': 1_000_000_000,
                         'min_eval_points': 2,
                         'units': {'x': 's', 'altitude_time': 'km', 'state_x_time': 'km',
                                   'state_y_time': 'km', 'state_z_time': 'km'},
                         'tolerances_quick': {'state_x_time': {'max_abs_km': 1e6}},
                         'tolerances_full': {'state_x_time': {'max_abs_km': 1e6}},
                         'spacecraft': {'bus_mass_kg': 29.0, 'bus_dims_m': [0.205, 0.37, 0.08],
                                        'panel_dims_m': [0.01, 0.0285, 0.0001]},
                         'atmosphere_truth': {'atmosphere_model': 'GRAM', 'gram_seed': 1001},
                         'calibration': {'enabled': False}, 'drag_enabled': False, 'EI_km': 300.0}
                if target == 'gmat' or degree == 0:
                    model['gravity_harmonics_gm_override_m3s2'] = 3.986004418e14
                self.write_manifest(case, model)
                gravity = {'path': f'/retained/repo/{field}', 'sha256': 'd'*64}
                kernel = {'path': '/retained/kernels/spk/planets/de430.bsp', 'sha256': 'e'*64}
                environment = {'SPACEAGORA_TELEMETRY_SOLVER_MODE': 'dp8',
                               'SPACEAGORA_TELEMETRY_DT_MAX_ORBIT': '10.0',
                               'SPACEAGORA_TELEMETRY_RELTOL_ORBIT': '1e-12',
                               'SPACEAGORA_TELEMETRY_ABSTOL_ORBIT': '1e-12',
                               'SPACEAGORA_TELEMETRY_RELTOL_ATM': '1e-12',
                               'SPACEAGORA_TELEMETRY_ABSTOL_ATM': '1e-12',
                               'SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH': 'spk/planets/de430.bsp'}
                cases[name] = {'inputs': [ref, gravity, kernel], 'kernels': [kernel], 'planetary_kernel': kernel,
                               'solver_environment': environment, 'model': model,
                               'solver_retcode': 'Success',
                               'manifest_sha256': export.digest(case / 'manifest.toml'),
                               'series_sha256': export.digest(case / 'series.arrow')}
                rows.append({'scenario': name, 'target': target, 'variant': 'committed',
                             'reference_path': ref['path'], 'solver_retcode': 'Success',
                             'rms_m': 1.0, 'max_m': 2.0, 'n_points': 3, 'arc_end_s': 2.0,
                             'rms_first10k_m': 1.0, 'gravity_degree': degree, 'gravity_order': degree,
                             'gravity_file': field,
                             'gm_override_m3s2': model.get('gravity_harmonics_gm_override_m3s2', 'NaN'),
                             'nbody_bodies': '+'.join(bodies)})
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

    @staticmethod
    def write_manifest(case_path, model):
        # A small writer for this fixture's primitive values and inline tables;
        # the exporter must parse actual TOML rather than trusting its checksum.
        def value(item):
            if isinstance(item, dict):
                return '{' + ', '.join(f'{json.dumps(k)} = {value(v)}' for k, v in item.items()) + '}'
            return json.dumps(item)
        text = '[[scenarios]]\n' + ''.join(f'{json.dumps(k)} = {value(v)}\n' for k, v in model.items())
        (case_path / 'manifest.toml').write_text(text)

    def edit_case(self, fn, *, name='earth_j50_tbfalse', rewrite_manifest=False):
        run = self.paths[0]
        info = json.loads((run / 'run_info.json').read_text())
        with (run / 'results.csv').open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        row = next(row for row in rows if row['scenario'] == name)
        case = info['cases'][name]
        fn(case, row)
        if rewrite_manifest:
            self.write_manifest(run / name, case['model'])
            case['manifest_sha256'] = export.digest(run / name / 'manifest.toml')
        with (run / 'results.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        info['results_sha256'] = export.digest(run / 'results.csv')
        self.save(run, info)

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

    def test_allowed_selection_and_warning_controls(self):
        self.mutate(lambda i: i.update(controls={'XVAL_SCENARIOS': ','.join(sorted(export.SCENARIOS)),
                                                'SPACEAGORA_WARN_NORMALIZE': 'false',
                                                'SPACEAGORA_WARN_DEPRECATED_CONFIG': 'false'}))
        self.assertEqual(len(export.validate_campaign(self.root)), 2)

    # These six independent counterexamples reproduce the published review.
    def test_missing_declared_gravity_input_is_rejected(self):  # C1
        self.edit_case(lambda case, row: case['inputs'].pop(1))
        self.rejected()
        destination = self.root / 'paper'
        with self.assertRaisesRegex(ValueError, 'gravity file input identity'):
            export.main(self.root, destination)
        self.assertFalse(destination.exists())

    def test_name_only_model_is_rejected(self):  # C2
        self.edit_case(lambda case, row: case.update(model={'name': row['scenario']}))
        self.rejected()

    def test_kernel_only_solver_environment_is_rejected(self):  # C3
        self.edit_case(lambda case, row: case.update(solver_environment={
            'SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH': 'spk/planets/de430.bsp'}))
        self.rejected()

    def test_model_gravity_path_must_match_manifest(self):  # C4
        self.edit_case(lambda case, row: case['model'].update(gravity_harmonics_file='other.csv'))
        with self.assertRaisesRegex(ValueError, 'Manifest and recorded model differ'):
            export.validate_campaign(self.root)

    def test_identity_for_a_different_gravity_file_is_rejected(self):  # C5
        self.edit_case(lambda case, row: case['inputs'][1].update(
            path='/retained/repo/data/Gravity_harmonics_data/other.csv'))
        with self.assertRaisesRegex(ValueError, 'gravity file input identity'):
            export.validate_campaign(self.root)

    def test_result_degree_must_match_model_and_manifest(self):  # C6
        self.edit_case(lambda case, row: case['model'].update(gravity_harmonics_degree=2,
                                                            gravity_harmonics_order=2),
                       rewrite_manifest=True)
        with self.assertRaisesRegex(ValueError, 'Result and model gravity degree/order differ'):
            export.validate_campaign(self.root)

    def test_each_substantive_model_field_is_required(self):
        original = json.loads((self.paths[0] / 'run_info.json').read_text())
        for key in ('name', 'planet', 'gravity_model', 'gravity_harmonics_degree',
                    'gravity_harmonics_order', 'gravity_harmonics_file', 'nbody_bodies',
                    'orbit_altitude_mode'):
            with self.subTest(key=key):
                info = copy.deepcopy(original)
                del info['cases']['earth_j50_tbfalse']['model'][key]
                self.save(self.paths[0], info)
                self.rejected()

    def test_base_model_fields_cannot_be_stripped_from_both_records(self):
        run, name = self.paths[0], 'earth_j50_tbfalse'
        original = json.loads((run / 'run_info.json').read_text())
        for key in ('kind', 'comparison_mode', 'events', 'telemetry', 'telemetry_columns',
                    'max_points_quick', 'max_points_full', 'min_eval_points', 'units',
                    'tolerances_quick', 'tolerances_full', 'initial_time', 'spacecraft',
                    'atmosphere_truth', 'calibration', 'drag_enabled', 'EI_km'):
            with self.subTest(key=key):
                info = copy.deepcopy(original)
                case = info['cases'][name]
                del case['model'][key]
                self.write_manifest(run / name, case['model'])
                case['manifest_sha256'] = export.digest(run / name / 'manifest.toml')
                self.save(run, info)
                self.rejected()

    def test_each_substantive_solver_field_is_required(self):
        original = json.loads((self.paths[0] / 'run_info.json').read_text())
        for key in ('SPACEAGORA_TELEMETRY_SOLVER_MODE', 'SPACEAGORA_TELEMETRY_DT_MAX_ORBIT',
                    'SPACEAGORA_TELEMETRY_RELTOL_ORBIT', 'SPACEAGORA_TELEMETRY_ABSTOL_ORBIT',
                    'SPACEAGORA_TELEMETRY_RELTOL_ATM', 'SPACEAGORA_TELEMETRY_ABSTOL_ATM',
                    'SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH'):
            with self.subTest(key=key):
                info = copy.deepcopy(original)
                del info['cases']['earth_j50_tbfalse']['solver_environment'][key]
                self.save(self.paths[0], info)
                self.rejected()

    def test_invalid_solver_settings_are_rejected(self):
        original = json.loads((self.paths[0] / 'run_info.json').read_text())
        for value in ('', 'NaN', 'Inf', '-1', '0', 'not-a-number', 10.0):
            with self.subTest(value=value):
                info = copy.deepcopy(original)
                info['cases']['earth_j50_tbfalse']['solver_environment']['SPACEAGORA_TELEMETRY_DT_MAX_ORBIT'] = value
                self.save(self.paths[0], info)
                self.rejected()

    def test_solver_kernel_selector_must_match_identity(self):
        self.edit_case(lambda case, row: case['solver_environment'].update(
            SPACEAGORA_SPICE_PLANETARY_KERNEL_RELPATH='spk/planets/other.bsp'))
        self.rejected()

    def test_relative_and_absolute_gravity_paths(self):
        # The original checkout and its inputs need not be available at export.
        self.assertEqual(len(export.validate_campaign(self.root)), 2)
        def absolute(case, row):
            field = case['inputs'][1]['path']
            case['model']['gravity_harmonics_file'] = field
            row['gravity_file'] = field
        self.edit_case(absolute, rewrite_manifest=True)
        self.assertEqual(len(export.validate_campaign(self.root)), 2)
        def normalized(case, row):
            field = 'data/./Gravity_harmonics_data/earth.csv'
            case['model']['gravity_harmonics_file'] = field
            row['gravity_file'] = field
        self.edit_case(normalized, rewrite_manifest=True)
        self.assertEqual(len(export.validate_campaign(self.root)), 2)

    def test_ambiguous_relative_gravity_identity_is_rejected(self):
        self.edit_case(lambda case, row: case['inputs'].append(
            {'path': '/another/repo/data/Gravity_harmonics_data/earth.csv', 'sha256': 'f'*64}))
        self.rejected()

    def test_absolute_path_must_match_whole_identity(self):
        def change(case, row):
            field = '/different/repo/data/Gravity_harmonics_data/earth.csv'
            case['model']['gravity_harmonics_file'] = field
            row['gravity_file'] = field
        self.edit_case(change, rewrite_manifest=True)
        self.rejected()

    def test_relative_parent_traversal_is_rejected(self):
        def change(case, row):
            field = 'data/../Gravity_harmonics_data/earth.csv'
            case['model']['gravity_harmonics_file'] = field
            row['gravity_file'] = field
        self.edit_case(change, rewrite_manifest=True)
        self.rejected()

    def test_point_mass_without_declared_file(self):
        def point_mass(case, row):
            case['model']['gravity_harmonics_file'] = ''
            del case['model']['gravity_harmonics_gm_override_m3s2']
            case['inputs'].pop(1)
            row.update(gravity_file='', gm_override_m3s2='NaN')
        self.edit_case(point_mass, name='earth_j0_tbfalse', rewrite_manifest=True)
        self.assertEqual(len(export.validate_campaign(self.root)), 2)

    def test_nonzero_degree_requires_declared_file(self):
        def change(case, row):
            case['model']['gravity_harmonics_file'] = ''
            del case['model']['gravity_harmonics_gm_override_m3s2']
            case['inputs'].pop(1)
            row.update(gravity_file='', gm_override_m3s2='NaN')
        self.edit_case(change, rewrite_manifest=True)
        with self.assertRaisesRegex(ValueError, 'Model requires a declared gravity file'):
            export.validate_campaign(self.root)

    def test_point_mass_gm_override_requires_declared_file(self):
        def change(case, row):
            case['model']['gravity_harmonics_file'] = ''
            case['inputs'].pop(1)
            row['gravity_file'] = ''
        self.edit_case(change, name='earth_j0_tbfalse', rewrite_manifest=True)
        with self.assertRaisesRegex(ValueError, 'Model requires a declared gravity file'):
            export.validate_campaign(self.root)

    def test_point_mass_with_declared_file_still_needs_identity(self):
        self.edit_case(lambda case, row: case['inputs'].pop(1), name='earth_j0_tbfalse')
        self.rejected()

    def test_result_configuration_fields_match_model(self):
        run = self.paths[0]
        original = (run / 'results.csv').read_bytes()
        for column, value in (('gravity_order', '0'), ('gravity_file', 'other.csv'),
                              ('gm_override_m3s2', '1.0'), ('nbody_bodies', 'sun')):
            with self.subTest(column=column):
                (run / 'results.csv').write_bytes(original)
                self.edit_case(lambda case, row: row.update({column: value}))
                self.rejected()

    def test_absent_gm_cannot_have_finite_result_override(self):
        def change(case, row):
            del case['model']['gravity_harmonics_gm_override_m3s2']
        self.edit_case(change, rewrite_manifest=True)
        self.rejected()

    def test_checksummed_manifest_must_be_valid_and_match_case(self):
        path = self.paths[0] / 'earth_j50_tbfalse' / 'manifest.toml'
        for content in ('this is not TOML', 'fixture = true\n', 'scenarios = []\n'):
            with self.subTest(content=content):
                path.write_text(content)
                self.mutate(lambda i: i['cases']['earth_j50_tbfalse'].update(manifest_sha256=export.digest(path)))
                self.rejected()

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
            self.assertEqual(row['gravity_file'], f"data/Gravity_harmonics_data/{row['body']}.csv")
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
