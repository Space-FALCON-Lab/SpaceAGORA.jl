"""Primary-export identity checks; no pandas, pyarrow, telemetry, or campaign."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "xval_fullarc_export.py"
spec = importlib.util.spec_from_file_location("fullarc_export", SCRIPT)
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)


def write_toml(path, data):
    lines = []
    def table(values, parents=()):
        for key, value in values.items():
            if not isinstance(value, dict):
                if isinstance(value, list) and value and isinstance(value[0], dict):
                    encoded = "[" + ", ".join("{" + ", ".join(f"{json.dumps(k)} = {json.dumps(v)}"
                                            for k, v in item.items()) + "}" for item in value) + "]"
                else:
                    encoded = json.dumps(value)
                lines.append(f"{json.dumps(key)} = {encoded}")
        for key, value in values.items():
            if isinstance(value, dict):
                parts = (*parents, key)
                lines.append("[" + ".".join(json.dumps(p) for p in parts) + "]")
                table(value, parts)
    table(data)
    path.write_text("\n".join(lines) + "\n")


class PrimaryProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name) / "gmat_committed"
        self.run.mkdir()
        self.scenarios = [f"{b}_{g}_{tb}" for b in export.BODY_ORDER for g in export.GRAV_ORDER
                          for tb in ("tbfalse", "tbtrue")]
        self.rows = [dict(scenario=s, target="gmat", variant="committed",
                          reference_sha256="a" * 64, gravity_file_sha256="b" * 64)
                     for s in self.scenarios]
        self.info = dict(provenance_schema=1, status="complete", target="gmat", variant="committed",
                         primary_eligible=True, input_overrides={}, src_test_scripts_dirty=False,
                         commit="c" * 40, scenarios=self.scenarios.copy(),
                         model_inputs=dict(planetary_kernel="fixture.bsp", planetary_kernel_sha256="d" * 64,
                                           solver_environment={"SPACEAGORA_TELEMETRY_SOLVER_MODE": "dp8"}),
                         spice_kernels=[dict(path="fixture.bsp", kind="SPK", sha256="d" * 64)],
                         artifacts_sha256={})
        for scenario in self.scenarios:
            (self.run / scenario).mkdir()
            for name in ("manifest.toml", "series.arrow"):
                path = self.run / scenario / name
                path.write_text("synthetic " + scenario)
                self.info["artifacts_sha256"][f"{scenario}/{name}"] = export.file_digest(path)
        self.save_rows()
        self.save_info()

    def save_rows(self):
        path = self.run / "results.csv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(self.rows[0]))
            writer.writeheader()
            writer.writerows(self.rows)
        self.info["artifacts_sha256"]["results.csv"] = export.file_digest(path)

    def save_info(self):
        write_toml(self.run / "run_info.toml", self.info)

    def validate(self):
        return export.validate_primary_run(self.run, "gmat")

    def test_complete_committed_identity_is_accepted(self):
        result = self.validate()
        self.assertEqual(result["commit"], "c" * 40)
        self.assertEqual(result["run_info_sha256"], export.file_digest(self.run / "run_info.toml"))

    def test_legacy_metadata_cannot_bless_primary(self):
        del self.info["provenance_schema"]
        self.save_info()
        with self.assertRaisesRegex(ValueError, "legacy provenance"):
            self.validate()
        (self.run / "run_info.toml").unlink()
        with self.assertRaisesRegex(ValueError, "missing versioned"):
            self.validate()

    def test_sensitive_or_unknown_inputs_are_not_primary(self):
        for key in ("XVAL_FRAME_TABLE_DIR", "XVAL_PLANETARY_KERNEL", "SPACEAGORA_SPICE_PCK_OVERRIDES",
                    "XVAL_BASE", "XVAL_FUTURE_MODEL_OVERRIDE", "SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT",
                    "SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT", "SPACEAGORA_GMAT_PARITY_SOLVER"):
            with self.subTest(key=key):
                self.info["input_overrides"] = {key: "diagnostic"}
                self.save_info()
                with self.assertRaisesRegex(ValueError, "overrides"):
                    self.validate()

    def test_renamed_sensitivity_directory_is_rejected(self):
        self.info["variant"] = "frame_diagnostic"
        self.save_info()
        with self.assertRaisesRegex(ValueError, "identity"):
            self.validate()

    def test_dirty_or_unidentified_sources_are_rejected(self):
        self.info["src_test_scripts_dirty"] = True
        self.save_info()
        with self.assertRaisesRegex(ValueError, "not clean"):
            self.validate()
        self.info["src_test_scripts_dirty"] = False
        self.info["commit"] = "unknown"
        self.save_info()
        with self.assertRaisesRegex(ValueError, "source commit"):
            self.validate()

    def test_interrupted_rerun_is_rejected_before_any_export(self):
        self.info["status"] = "running"
        self.save_info()
        paper = Path(self.tmp.name) / "paper"
        with self.assertRaisesRegex(ValueError, "not complete"):
            export.main(self.tmp.name, paper)
        self.assertFalse(paper.exists())

    def test_modified_results_or_series_are_rejected(self):
        for relative in ("results.csv", self.scenarios[0] + "/series.arrow"):
            with self.subTest(artifact=relative):
                path = self.run / relative
                original = path.read_bytes()
                path.write_bytes(original + b"changed")
                with self.assertRaisesRegex(ValueError, "artifact differs"):
                    self.validate()
                path.write_bytes(original)

    def test_subset_or_duplicate_selection_is_not_full_matrix(self):
        for selected in (self.scenarios[:-1], self.scenarios[:-1] + [self.scenarios[0]]):
            self.info["scenarios"] = selected
            self.save_info()
            with self.assertRaisesRegex(ValueError, "24 scenarios exactly once"):
                self.validate()

    def test_rehashed_duplicate_rows_are_still_rejected(self):
        self.rows[-1] = self.rows[0].copy()
        self.save_rows()
        self.save_info()
        with self.assertRaisesRegex(ValueError, "complete unique"):
            self.validate()

    def test_row_target_and_input_identity_are_checked(self):
        self.rows[0]["target"] = "stk"
        self.save_rows()
        self.save_info()
        with self.assertRaisesRegex(ValueError, "different target"):
            self.validate()
        self.rows[0]["target"] = "gmat"
        self.rows[0]["reference_sha256"] = ""
        self.save_rows()
        self.save_info()
        with self.assertRaisesRegex(ValueError, "reference or gravity"):
            self.validate()

    def test_loaded_kernel_identity_is_required(self):
        self.info["spice_kernels"][0]["sha256"] = "e" * 64
        self.save_info()
        with self.assertRaisesRegex(ValueError, "planetary-kernel mismatch"):
            self.validate()

    def test_missing_effective_inputs_or_artifact_map_fails_closed(self):
        inputs = self.info.pop("model_inputs")
        self.save_info()
        with self.assertRaisesRegex(ValueError, "effective kernel"):
            self.validate()
        self.info["model_inputs"] = inputs
        self.info["artifacts_sha256"].pop("results.csv")
        self.save_info()
        with self.assertRaisesRegex(ValueError, "artifact identity"):
            self.validate()


if __name__ == "__main__":
    unittest.main()
