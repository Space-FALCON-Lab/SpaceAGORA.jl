from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orchestrator.app.config import Settings
from orchestrator.app.models import (
    FilePaths,
    InitialConditions,
    OutputOptions,
    SimulationSpec,
    SimulationModels,
)
from orchestrator.app.services.capability_validator import CapabilityValidator
from orchestrator.app.models import RunMode


def _make_settings(root: Path) -> Settings:
    return Settings(
        repo_root=root,
        storage_dir=root / "data",
        runs_dir=root / "data" / "runs",
        artifacts_dir=root / "data" / "artifacts",
        worker_tmp_dir=root / "tmp",
        spice_path=root / "GRAM_Data" / "SPICE",
        gram_data_path=root / "GRAM_Data",
        gram_py_path=root / "data/GRAMSuite.jl/GRAM Suite 2.0",
        max_mission_time_sec=604800.0,
        sync_max_mission_time_sec=21600.0,
        max_sync_job_wall_sec=180.0,
        enable_web_retrieval=False,
        default_run_mode="async",
    )


def _create_spice_tree(base: Path) -> None:
    (base / "pck").mkdir(parents=True, exist_ok=True)
    (base / "lsk").mkdir(parents=True, exist_ok=True)
    (base / "spk" / "planets").mkdir(parents=True, exist_ok=True)
    (base / "pck" / "pck00011.tpc").write_text("x", encoding="utf-8")
    (base / "lsk" / "naif0012.tls").write_text("x", encoding="utf-8")
    (base / "spk" / "planets" / "de440s.bsp").write_text("x", encoding="utf-8")


class CapabilityValidatorTests(unittest.TestCase):
    def test_rejects_near_circular_orbit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _create_spice_tree(root / "GRAM_Data" / "SPICE")
            (root / "GRAM_Data").mkdir(parents=True, exist_ok=True)
            (root / "data/GRAMSuite.jl/GRAM Suite 2.0").mkdir(parents=True, exist_ok=True)
            settings = _make_settings(root)
            validator = CapabilityValidator(settings)

            spec = SimulationSpec(
                run_name="test",
                planet="Mars",
                mission_time_sec=1000.0,
                entry_interface_km=125.0,
                initial_conditions=InitialConditions(
                    ra_m=7000000.0,
                    rp_m=6999999.0,
                    i_deg=30.0,
                    argument_of_periapsis_deg=0.0,
                    raan_deg=0.0,
                    true_anomaly_deg=180.0,
                ),
                models=SimulationModels(density_model="gram"),
                outputs=OutputOptions(results=True, save_csv=True),
                file_paths=FilePaths(
                    spice=str(root / "GRAM_Data" / "SPICE"),
                    gram_data=str(root / "GRAM_Data"),
                    gram_py=str(root / "data/GRAMSuite.jl/GRAM Suite 2.0"),
                ),
            )

            _, report = validator.validate(spec, requested_mode=RunMode.ASYNC)
            self.assertFalse(report.supported)
            self.assertTrue(any("near-circular" in r for r in report.unsupported_reasons))

    def test_forces_csv_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _create_spice_tree(root / "GRAM_Data" / "SPICE")
            (root / "GRAM_Data").mkdir(parents=True, exist_ok=True)
            (root / "data/GRAMSuite.jl/GRAM Suite 2.0").mkdir(parents=True, exist_ok=True)
            settings = _make_settings(root)
            validator = CapabilityValidator(settings)

            spec = SimulationSpec(
                run_name="test",
                planet="Earth",
                mission_time_sec=1000.0,
                entry_interface_km=120.0,
                initial_conditions=InitialConditions(
                    ra_m=7500000.0,
                    rp_m=6800000.0,
                    i_deg=30.0,
                    argument_of_periapsis_deg=0.0,
                    raan_deg=0.0,
                    true_anomaly_deg=180.0,
                ),
                models=SimulationModels(density_model="none"),
                outputs=OutputOptions(results=False, save_csv=False),
                file_paths=FilePaths(
                    spice=str(root / "GRAM_Data" / "SPICE"),
                    gram_data=str(root / "GRAM_Data"),
                    gram_py=str(root / "data/GRAMSuite.jl/GRAM Suite 2.0"),
                ),
            )

            validated, report = validator.validate(spec, requested_mode=RunMode.SYNC)
            self.assertTrue(report.supported)
            self.assertTrue(validated.outputs.results)
            self.assertTrue(validated.outputs.save_csv)
            self.assertEqual(report.forced_overrides["outputs.results"], True)
            self.assertEqual(report.forced_overrides["outputs.save_csv"], True)


if __name__ == "__main__":
    unittest.main()
