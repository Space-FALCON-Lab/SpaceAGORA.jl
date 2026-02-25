from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, HttpUrl


class RunMode(str, Enum):
    SYNC = "sync"
    ASYNC = "async"


class RunState(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    UNSUPPORTED = "unsupported"


class SourceCitation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    label: str
    url: str
    retrieved_at: datetime | None = None
    notes: str | None = None


class MissionIntent(BaseModel):
    model_config = ConfigDict(extra="forbid")

    objective: str
    mission_name: str | None = None
    planet: str | None = None
    replicate_mission: bool = False
    requested_constraints: dict[str, Any] = Field(default_factory=dict)


class InitialConditions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ra_m: float
    rp_m: float
    i_deg: float
    omega_deg: float = Field(alias="argument_of_periapsis_deg")
    raan_deg: float
    true_anomaly_deg: float = 180.0


class SimulationModels(BaseModel):
    model_config = ConfigDict(extra="forbid")

    gravity_model: str = "inverse_squared"
    density_model: str = "gram"
    thermal_model: str = "maxwellian"
    aerodynamic_model: str = "mach_dependent"
    wind: bool = True
    srp: bool = False
    n_bodies: list[str] = Field(default_factory=list)


class IntegrationTolerances(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reltol_orbit: float = 1e-6
    abstol_orbit: float = 1e-7
    reltol_atmosphere: float = 1e-7
    abstol_atmosphere: float = 1e-9
    dt_max_orbit: float = 30.0
    dt_max_atmosphere: float = 1.0


class OutputOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    results: bool = True
    save_csv: bool = True
    generate_plots: bool = False
    include_log_bundle: bool = True


class FilePaths(BaseModel):
    model_config = ConfigDict(extra="forbid")

    spice: str
    gram_data: str
    gram_py: str


class SimulationSpec(BaseModel):
    model_config = ConfigDict(extra="forbid")

    spec_version: str = "1.0"
    run_name: str
    planet: str
    mission_name: str | None = None
    mission_time_sec: float
    entry_interface_km: float = 125.0
    initial_conditions: InitialConditions
    models: SimulationModels = Field(default_factory=SimulationModels)
    tolerances: IntegrationTolerances = Field(default_factory=IntegrationTolerances)
    outputs: OutputOptions = Field(default_factory=OutputOptions)
    file_paths: FilePaths
    keplerian: bool = False
    orientation_sim: bool = False
    num_steps_to_save: int = 1000


class CapabilityReport(BaseModel):
    model_config = ConfigDict(extra="forbid")

    supported: bool
    unsupported_reasons: list[str] = Field(default_factory=list)
    nearest_supported_spec: SimulationSpec | None = None
    forced_overrides: dict[str, Any] = Field(default_factory=dict)
    notes: list[str] = Field(default_factory=list)
    reroute_to_async: bool = False


class ArtifactRef(BaseModel):
    model_config = ConfigDict(extra="forbid")

    artifact_id: str
    kind: str
    path: str
    content_type: str
    size_bytes: int | None = None


class ResultSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    summary_text: str
    key_metrics: dict[str, Any] = Field(default_factory=dict)
    caveats: list[str] = Field(default_factory=list)
    citations: list[SourceCitation] = Field(default_factory=list)
    artifacts: list[ArtifactRef] = Field(default_factory=list)


class RunRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    status: RunState
    requested_mode: RunMode
    effective_mode: RunMode
    prompt: str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    error_message: str | None = None
    intent: MissionIntent | None = None
    spec: SimulationSpec | None = None
    capability_report: CapabilityReport | None = None
    result: ResultSummary | None = None


class CreateRunRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    prompt: str = Field(min_length=3)
    mode: RunMode | None = None
    mission_overrides: dict[str, Any] = Field(default_factory=dict)
    output_preferences: dict[str, Any] = Field(default_factory=dict)


class CreateRunResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    status: RunState
    accepted: bool
    requested_mode: RunMode
    effective_mode: RunMode
    capability_report: CapabilityReport


class RunStatusResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    status: RunState
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    error_message: str | None = None


class RunResultResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    status: RunState
    result: ResultSummary


class ErrorResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    detail: str


def now_utc() -> datetime:
    return datetime.now(timezone.utc)
