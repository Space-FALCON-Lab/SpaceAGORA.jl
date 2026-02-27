from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    repo_root: Path
    storage_dir: Path
    runs_dir: Path
    artifacts_dir: Path
    worker_tmp_dir: Path
    spice_path: Path
    gram_data_path: Path
    gram_py_path: Path
    max_mission_time_sec: float
    sync_max_mission_time_sec: float
    max_sync_job_wall_sec: float
    enable_web_retrieval: bool
    default_run_mode: str


def _bool_env(key: str, default: bool) -> bool:
    value = os.getenv(key)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _first_existing_path(*paths: Path) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]


def load_settings() -> Settings:
    # This file sits at orchestrator/app/config.py
    repo_root = Path(__file__).resolve().parents[2]
    storage_dir = Path(os.getenv("ORCH_STORAGE_DIR", str(repo_root / "orchestrator" / "data")))
    runs_dir = Path(os.getenv("ORCH_RUNS_DIR", str(storage_dir / "runs")))
    artifacts_dir = Path(os.getenv("ORCH_ARTIFACTS_DIR", str(storage_dir / "artifacts")))
    worker_tmp_dir = Path(os.getenv("ORCH_WORKER_TMP_DIR", "/tmp/spaceagora_orchestrator"))

    default_spice = _first_existing_path(repo_root / "GRAM Suite 2.0" / "SPICE", repo_root / "GRAM_Data" / "SPICE")
    default_gram_data = _first_existing_path(repo_root / "GRAM Suite 2.0", repo_root / "GRAM_Data")
    default_gram_root = _first_existing_path(repo_root / "GRAM Suite 2.0", repo_root / "GRAMpy")

    spice_path = Path(os.getenv("ORCH_SPICE_PATH", str(default_spice)))
    gram_data_path = Path(os.getenv("ORCH_GRAM_DATA_PATH", str(default_gram_data)))
    gram_root_env = os.getenv("ORCH_GRAM_ROOT_PATH")
    gram_py_path = Path(os.getenv("ORCH_GRAMPY_PATH", gram_root_env if gram_root_env is not None else str(default_gram_root)))

    max_mission_time_sec = float(os.getenv("ORCH_MAX_MISSION_TIME_SEC", "604800"))
    sync_max_mission_time_sec = float(os.getenv("ORCH_SYNC_MAX_MISSION_TIME_SEC", "21600"))
    max_sync_job_wall_sec = float(os.getenv("ORCH_MAX_SYNC_JOB_WALL_SEC", "180"))

    enable_web_retrieval = _bool_env("ORCH_ENABLE_WEB_RETRIEVAL", False)
    default_run_mode = os.getenv("ORCH_DEFAULT_RUN_MODE", "async")

    runs_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    worker_tmp_dir.mkdir(parents=True, exist_ok=True)

    return Settings(
        repo_root=repo_root,
        storage_dir=storage_dir,
        runs_dir=runs_dir,
        artifacts_dir=artifacts_dir,
        worker_tmp_dir=worker_tmp_dir,
        spice_path=spice_path,
        gram_data_path=gram_data_path,
        gram_py_path=gram_py_path,
        max_mission_time_sec=max_mission_time_sec,
        sync_max_mission_time_sec=sync_max_mission_time_sec,
        max_sync_job_wall_sec=max_sync_job_wall_sec,
        enable_web_retrieval=enable_web_retrieval,
        default_run_mode=default_run_mode,
    )
