from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orchestrator.app.adapters.artifacts import LocalArtifactStore
from orchestrator.app.adapters.queues import InMemoryJobQueue
from orchestrator.app.adapters.stores import FileRunStore
from orchestrator.app.config import Settings
from orchestrator.app.models import CreateRunRequest, RunMode, RunState
from orchestrator.app.services.capability_validator import CapabilityValidator
from orchestrator.app.services.intent_extractor import RuleBasedIntentExtractor
from orchestrator.app.services.mission_catalog import CuratedMissionCatalog
from orchestrator.app.services.mission_retriever import MissionRetriever
from orchestrator.app.services.run_orchestrator import RunOrchestrator
from orchestrator.app.services.spec_compiler import SimulationSpecCompiler
from orchestrator.app.worker.julia_runner import ExecutionResult


class FakeRunner:
    def __init__(self, root: Path) -> None:
        self._root = root

    def run(self, job_id: str, spec, wall_timeout_sec: float) -> ExecutionResult:
        run_dir = self._root / "worker" / job_id
        run_dir.mkdir(parents=True, exist_ok=True)

        spec_snapshot = run_dir / "simulation_spec.json"
        spec_snapshot.write_text(spec.model_dump_json(indent=2), encoding="utf-8")

        script_path = run_dir / "run_job.jl"
        script_path.write_text("println(\"fake\")\n", encoding="utf-8")

        stdout_path = run_dir / "stdout.log"
        stdout_path.write_text("RUN_OK\n", encoding="utf-8")

        stderr_path = run_dir / "stderr.log"
        stderr_path.write_text("", encoding="utf-8")

        csv_path = run_dir / "simulation_results.csv"
        csv_path.write_text(
            "time,sc1_pos_1,sc1_pos_2,sc1_pos_3,sc1_vel_1,sc1_vel_2,sc1_vel_3\n"
            "0.0,1,2,3,4,5,6\n"
            "10.0,7,8,9,10,11,12\n",
            encoding="utf-8",
        )

        return ExecutionResult(
            success=True,
            run_dir=run_dir,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            spec_snapshot_path=spec_snapshot,
            script_path=script_path,
            result_csv_path=csv_path,
        )


def _make_settings(root: Path) -> Settings:
    runs_dir = root / "data" / "runs"
    artifacts_dir = root / "data" / "artifacts"
    worker_tmp = root / "tmp"
    runs_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    worker_tmp.mkdir(parents=True, exist_ok=True)
    return Settings(
        repo_root=root,
        storage_dir=root / "data",
        runs_dir=runs_dir,
        artifacts_dir=artifacts_dir,
        worker_tmp_dir=worker_tmp,
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


class RunOrchestratorTests(unittest.TestCase):
    def test_unsupported_planet_request(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _create_spice_tree(root / "GRAM_Data" / "SPICE")
            (root / "data/GRAMSuite.jl/GRAM Suite 2.0").mkdir(parents=True, exist_ok=True)

            settings = _make_settings(root)
            orchestrator = RunOrchestrator(
                settings=settings,
                store=FileRunStore(settings.runs_dir),
                queue=InMemoryJobQueue(),
                artifact_store=LocalArtifactStore(settings.artifacts_dir),
                intent_extractor=RuleBasedIntentExtractor(),
                mission_retriever=MissionRetriever(CuratedMissionCatalog(), enable_web_fallback=False),
                spec_compiler=SimulationSpecCompiler(settings),
                validator=CapabilityValidator(settings),
                simulation_runner=FakeRunner(root),
            )

            request = CreateRunRequest(prompt="simulate a satellite around Venus", mode=RunMode.ASYNC)
            record = orchestrator.create_run(request)
            self.assertEqual(record.status, RunState.UNSUPPORTED)

    def test_sync_run_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _create_spice_tree(root / "GRAM_Data" / "SPICE")
            (root / "data/GRAMSuite.jl/GRAM Suite 2.0").mkdir(parents=True, exist_ok=True)

            settings = _make_settings(root)
            orchestrator = RunOrchestrator(
                settings=settings,
                store=FileRunStore(settings.runs_dir),
                queue=InMemoryJobQueue(),
                artifact_store=LocalArtifactStore(settings.artifacts_dir),
                intent_extractor=RuleBasedIntentExtractor(),
                mission_retriever=MissionRetriever(CuratedMissionCatalog(), enable_web_fallback=False),
                spec_compiler=SimulationSpecCompiler(settings),
                validator=CapabilityValidator(settings),
                simulation_runner=FakeRunner(root),
            )

            request = CreateRunRequest(prompt="replicate Mars Odyssey for 0.1 hours", mode=RunMode.SYNC)
            record = orchestrator.create_run(request)
            self.assertEqual(record.status, RunState.SUCCEEDED)
            self.assertIsNotNone(record.result)
            self.assertGreaterEqual(len(record.result.artifacts), 1)


if __name__ == "__main__":
    unittest.main()
