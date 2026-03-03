from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orchestrator.app.config import Settings
from orchestrator.app.models import MissionIntent
from orchestrator.app.services.mission_catalog import CuratedMissionCatalog
from orchestrator.app.services.mission_retriever import MissionRetriever
from orchestrator.app.services.spec_compiler import SimulationSpecCompiler


def _make_settings(root: Path) -> Settings:
    return Settings(
        repo_root=root,
        storage_dir=root / "data",
        runs_dir=root / "data" / "runs",
        artifacts_dir=root / "data" / "artifacts",
        worker_tmp_dir=root / "tmp",
        spice_path=root / "GRAM_Data" / "SPICE",
        gram_data_path=root / "GRAM_Data",
        gram_py_path=root / "GRAMSuite.jl/GRAM Suite 2.0",
        max_mission_time_sec=604800.0,
        sync_max_mission_time_sec=21600.0,
        max_sync_job_wall_sec=180.0,
        enable_web_retrieval=False,
        default_run_mode="async",
    )


class SimulationSpecCompilerTests(unittest.TestCase):
    def test_compiles_odyssey_template_with_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            settings = _make_settings(root)

            intent = MissionIntent(
                objective="simulate_mission",
                mission_name="Mars Odyssey",
                planet="Mars",
                replicate_mission=True,
            )
            retriever = MissionRetriever(CuratedMissionCatalog(), enable_web_fallback=False)
            retrieved = retriever.retrieve(intent)

            compiler = SimulationSpecCompiler(settings)
            spec, citations, notes = compiler.compile(
                job_id="abc123",
                intent=intent,
                retrieved=retrieved,
                mission_overrides={"mission_time_sec": 3600.0, "models": {"srp": True}},
                output_preferences={"generate_plots": False},
            )

            self.assertEqual(spec.planet, "Mars")
            self.assertEqual(spec.mission_name, "Mars Odyssey")
            self.assertEqual(spec.mission_time_sec, 3600.0)
            self.assertTrue(spec.models.srp)
            self.assertTrue(spec.outputs.results)
            self.assertTrue(spec.outputs.save_csv)
            self.assertGreaterEqual(len(citations), 1)
            self.assertIn("compiled_to_simulation_spec", notes)


if __name__ == "__main__":
    unittest.main()
