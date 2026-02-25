from __future__ import annotations

import uuid
from typing import Optional

from orchestrator.app.adapters.artifacts import ArtifactStore
from orchestrator.app.adapters.queues import JobQueue, RunJob
from orchestrator.app.adapters.stores import RunStore
from orchestrator.app.config import Settings
from orchestrator.app.errors import NotFoundError
from orchestrator.app.models import (
    ArtifactRef,
    CreateRunRequest,
    RunMode,
    RunRecord,
    RunState,
    now_utc,
)
from orchestrator.app.services.capability_validator import CapabilityValidator
from orchestrator.app.services.intent_extractor import IntentExtractor
from orchestrator.app.services.mission_retriever import MissionRetriever
from orchestrator.app.services.result_summarizer import summarize_result
from orchestrator.app.services.spec_compiler import SimulationSpecCompiler
from orchestrator.app.worker.julia_runner import JuliaSimulationRunner


class RunOrchestrator:
    def __init__(
        self,
        settings: Settings,
        store: RunStore,
        queue: JobQueue,
        artifact_store: ArtifactStore,
        intent_extractor: IntentExtractor,
        mission_retriever: MissionRetriever,
        spec_compiler: SimulationSpecCompiler,
        validator: CapabilityValidator,
        simulation_runner: JuliaSimulationRunner,
    ) -> None:
        self._settings = settings
        self._store = store
        self._queue = queue
        self._artifact_store = artifact_store
        self._intent_extractor = intent_extractor
        self._mission_retriever = mission_retriever
        self._spec_compiler = spec_compiler
        self._validator = validator
        self._simulation_runner = simulation_runner

    def create_run(self, request: CreateRunRequest) -> RunRecord:
        job_id = uuid.uuid4().hex
        requested_mode = request.mode or RunMode(self._settings.default_run_mode)

        intent = self._intent_extractor.extract(request.prompt)
        retrieved = self._mission_retriever.retrieve(intent)
        spec, citations, compile_notes = self._spec_compiler.compile(
            job_id=job_id,
            intent=intent,
            retrieved=retrieved,
            mission_overrides=request.mission_overrides,
            output_preferences=request.output_preferences,
        )
        spec, capability_report = self._validator.validate(spec, requested_mode=requested_mode)
        capability_report.notes.extend(compile_notes)

        effective_mode = RunMode.ASYNC if capability_report.reroute_to_async else requested_mode
        status = RunState.UNSUPPORTED if not capability_report.supported else RunState.QUEUED

        record = RunRecord(
            job_id=job_id,
            status=status,
            requested_mode=requested_mode,
            effective_mode=effective_mode,
            prompt=request.prompt,
            created_at=now_utc(),
            updated_at=now_utc(),
            intent=intent,
            spec=spec,
            capability_report=capability_report,
            citations=citations,
        )
        self._store.put(record)

        if not capability_report.supported:
            return record

        if effective_mode == RunMode.ASYNC:
            self._queue.enqueue(RunJob(job_id=job_id))
            return self._store.get(job_id)

        self.execute_job(job_id)
        return self._store.get(job_id)

    def get_run(self, job_id: str) -> RunRecord:
        return self._store.get(job_id)

    def execute_job(self, job_id: str) -> RunRecord:
        record = self._store.get(job_id)
        if record.spec is None:
            raise RuntimeError(f"Run {job_id} has no SimulationSpec")

        started_at = now_utc()
        record = self._store.update(
            job_id,
            lambda current: current.model_copy(
                update={
                    "status": RunState.RUNNING,
                    "started_at": started_at,
                    "updated_at": now_utc(),
                }
            ),
        )

        timeout = self._settings.max_sync_job_wall_sec
        if record.effective_mode == RunMode.ASYNC:
            timeout = max(timeout, 3600.0)

        execution = self._simulation_runner.run(job_id=job_id, spec=record.spec, wall_timeout_sec=timeout)

        artifacts: list[ArtifactRef] = []
        artifacts.append(self._artifact_store.save_file(job_id, execution.spec_snapshot_path, "spec", "config"))
        artifacts.append(self._artifact_store.save_file(job_id, execution.script_path, "julia_script", "script"))
        artifacts.append(self._artifact_store.save_file(job_id, execution.stdout_path, "stdout", "log"))
        artifacts.append(self._artifact_store.save_file(job_id, execution.stderr_path, "stderr", "log"))

        if execution.result_csv_path and execution.result_csv_path.exists():
            artifacts.append(
                self._artifact_store.save_file(
                    job_id,
                    execution.result_csv_path,
                    "trajectory",
                    "trajectory_csv",
                )
            )

        finished_at = now_utc()
        if execution.success:
            result = summarize_result(
                spec=record.spec,
                result_csv=execution.result_csv_path,
                citations=record.citations,
                artifacts=artifacts,
            )
            updated = self._store.update(
                job_id,
                lambda current: current.model_copy(
                    update={
                        "status": RunState.SUCCEEDED,
                        "finished_at": finished_at,
                        "updated_at": now_utc(),
                        "result": result,
                        "error_message": None,
                    }
                ),
            )
            return updated

        result = summarize_result(
            spec=record.spec,
            result_csv=execution.result_csv_path,
            citations=record.citations,
            artifacts=artifacts,
            extra_caveats=[execution.error_message or "simulation failed"],
        )
        updated = self._store.update(
            job_id,
            lambda current: current.model_copy(
                update={
                    "status": RunState.FAILED,
                    "finished_at": finished_at,
                    "updated_at": now_utc(),
                    "result": result,
                    "error_message": execution.error_message,
                }
            ),
        )
        return updated

    def execute_next_queued_job(self) -> Optional[RunRecord]:
        job = self._queue.dequeue(timeout_sec=0.5)
        if job is None:
            return None
        return self.execute_job(job.job_id)
