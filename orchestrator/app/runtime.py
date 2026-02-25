from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

from orchestrator.app.adapters.artifacts import ArtifactStore, LocalArtifactStore
from orchestrator.app.adapters.queues import InMemoryJobQueue, JobQueue
from orchestrator.app.adapters.stores import FileRunStore, RunStore
from orchestrator.app.config import Settings, load_settings
from orchestrator.app.services.capability_validator import CapabilityValidator
from orchestrator.app.services.intent_extractor import RuleBasedIntentExtractor
from orchestrator.app.services.mission_catalog import CuratedMissionCatalog
from orchestrator.app.services.mission_retriever import MissionRetriever
from orchestrator.app.services.run_orchestrator import RunOrchestrator
from orchestrator.app.services.spec_compiler import SimulationSpecCompiler
from orchestrator.app.services.worker_loop import BackgroundWorker
from orchestrator.app.worker.julia_runner import JuliaSimulationRunner


@dataclass
class AppRuntime:
    settings: Settings
    store: RunStore
    queue: JobQueue
    artifact_store: ArtifactStore
    orchestrator: RunOrchestrator
    worker: BackgroundWorker


_runtime_singleton: Optional[AppRuntime] = None


def build_runtime() -> AppRuntime:
    settings = load_settings()

    store = FileRunStore(settings.runs_dir)
    queue = InMemoryJobQueue()
    artifact_store = LocalArtifactStore(settings.artifacts_dir)

    intent_extractor = RuleBasedIntentExtractor()
    mission_catalog = CuratedMissionCatalog()
    mission_retriever = MissionRetriever(
        catalog=mission_catalog,
        enable_web_fallback=settings.enable_web_retrieval,
    )
    spec_compiler = SimulationSpecCompiler(settings)
    validator = CapabilityValidator(settings)
    simulation_runner = JuliaSimulationRunner(settings)

    orchestrator = RunOrchestrator(
        settings=settings,
        store=store,
        queue=queue,
        artifact_store=artifact_store,
        intent_extractor=intent_extractor,
        mission_retriever=mission_retriever,
        spec_compiler=spec_compiler,
        validator=validator,
        simulation_runner=simulation_runner,
    )
    worker = BackgroundWorker(orchestrator)

    return AppRuntime(
        settings=settings,
        store=store,
        queue=queue,
        artifact_store=artifact_store,
        orchestrator=orchestrator,
        worker=worker,
    )


def get_runtime() -> AppRuntime:
    global _runtime_singleton
    if _runtime_singleton is None:
        _runtime_singleton = build_runtime()
    return _runtime_singleton
