from __future__ import annotations

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse

from orchestrator.app.errors import NotFoundError
from orchestrator.app.models import (
    CreateRunRequest,
    CreateRunResponse,
    RunResultResponse,
    RunState,
    RunStatusResponse,
)
from orchestrator.app.runtime import get_runtime


app = FastAPI(title="SpaceAGORA LLM Orchestrator", version="0.1.0")


@app.on_event("startup")
def _startup() -> None:
    runtime = get_runtime()
    runtime.worker.start()


@app.on_event("shutdown")
def _shutdown() -> None:
    runtime = get_runtime()
    runtime.worker.stop()


@app.get("/healthz")
def healthz() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/v1/runs", response_model=CreateRunResponse)
def create_run(request: CreateRunRequest) -> CreateRunResponse:
    runtime = get_runtime()
    record = runtime.orchestrator.create_run(request)
    accepted = record.status != RunState.UNSUPPORTED
    return CreateRunResponse(
        job_id=record.job_id,
        status=record.status,
        accepted=accepted,
        requested_mode=record.requested_mode,
        effective_mode=record.effective_mode,
        capability_report=record.capability_report,
    )


@app.get("/v1/runs/{job_id}", response_model=RunStatusResponse)
def get_run_status(job_id: str) -> RunStatusResponse:
    runtime = get_runtime()
    try:
        record = runtime.orchestrator.get_run(job_id)
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    return RunStatusResponse(
        job_id=record.job_id,
        status=record.status,
        created_at=record.created_at,
        updated_at=record.updated_at,
        started_at=record.started_at,
        finished_at=record.finished_at,
        error_message=record.error_message,
    )


@app.get("/v1/runs/{job_id}/result", response_model=RunResultResponse)
def get_run_result(job_id: str) -> RunResultResponse:
    runtime = get_runtime()
    try:
        record = runtime.orchestrator.get_run(job_id)
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    if record.result is None:
        raise HTTPException(status_code=409, detail=f"Run '{job_id}' has no result yet")

    return RunResultResponse(job_id=record.job_id, status=record.status, result=record.result)


@app.get("/v1/runs/{job_id}/artifacts/{artifact_id}")
def get_artifact(job_id: str, artifact_id: str) -> FileResponse:
    runtime = get_runtime()

    try:
        record = runtime.orchestrator.get_run(job_id)
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    if record.result is None:
        raise HTTPException(status_code=409, detail=f"Run '{job_id}' has no artifacts yet")

    artifact = next((a for a in record.result.artifacts if a.artifact_id == artifact_id), None)
    if artifact is None:
        raise HTTPException(status_code=404, detail=f"Artifact '{artifact_id}' not found")

    path = runtime.artifact_store.resolve_path(job_id, artifact_id)
    return FileResponse(path=str(path), media_type=artifact.content_type, filename=path.name)
