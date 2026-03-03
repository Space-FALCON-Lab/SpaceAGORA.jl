# SpaceAGORA Cloud Orchestrator

Standalone service layer for prompt-driven mission simulation that keeps `SpaceAGORA.jl` solver code read-only.

## What This Implements
- FastAPI service with required endpoints:
  - `POST /v1/runs`
  - `GET /v1/runs/{job_id}`
  - `GET /v1/runs/{job_id}/result`
  - `GET /v1/runs/{job_id}/artifacts/{artifact_id}`
- Strict typed contracts (`MissionIntent`, `SimulationSpec`, `CapabilityReport`, `RunRecord`, `ResultSummary`)
- Hybrid mission retrieval (curated catalog first, optional web fallback)
- Capability validator enforcing known runtime constraints
- Async queue + background worker
- Julia worker bridge that executes `run_analysis(args::SimulationConfiguration)` without modifying `src/`
- Artifact packaging (spec snapshot, generated Julia script, stdout/stderr, trajectory CSV)

## Directory Layout
- `app/main.py`: API
- `app/models.py`: request/response + internal schemas
- `app/services/`: intent extraction, mission retrieval, spec compilation, validation, orchestration
- `app/worker/julia_runner.py`: Julia subprocess execution bridge
- `app/adapters/`: queue/store/artifact interfaces and local implementations
- `scripts/export_schemas.py`: export JSON schema files
- `tests/`: unit tests
- `infra/terraform/`: AWS deployment scaffolding

## Quick Start
1. Install dependencies:
   ```bash
   pip install -r orchestrator/requirements.txt
   ```
2. Start API:
   ```bash
   uvicorn orchestrator.app.main:app --host 0.0.0.0 --port 8080
   ```
3. Submit a run:
   ```bash
   curl -X POST http://localhost:8080/v1/runs \
     -H 'Content-Type: application/json' \
     -d '{"prompt":"replicate Mars Odyssey", "mode":"async"}'
   ```

## Environment Variables
- `ORCH_STORAGE_DIR`: run metadata and artifact root
- `ORCH_SPICE_PATH`: SPICE kernel root (default: `<repo>/GRAMSuite.jl/GRAM Suite 2.0/SPICE`, fallback `<repo>/GRAM_Data/SPICE`)
- `ORCH_GRAM_DATA_PATH`: GRAM data root (default: `<repo>/GRAMSuite.jl/GRAM Suite 2.0`, fallback `<repo>/GRAM_Data`)
- `ORCH_GRAM_ROOT_PATH`: GRAM Suite root for Julia wrapper (legacy alias: `ORCH_GRAMPY_PATH`)
- `ORCH_ENABLE_WEB_RETRIEVAL`: `true/false` (default `false`)
- `ORCH_MAX_MISSION_TIME_SEC`: hard mission-time cap
- `ORCH_SYNC_MAX_MISSION_TIME_SEC`: sync-mode mission-time cap before async reroute
- `ORCH_MAX_SYNC_JOB_WALL_SEC`: sync execution timeout

## Current Notes
- Validator force-enables `outputs.results=true` and `outputs.save_csv=true` due known runner behavior.
- Validator blocks near-circular initial orbits to avoid callback/OE instability.
- Curated Mars catalog includes Mars Odyssey template for deterministic launch behavior.
