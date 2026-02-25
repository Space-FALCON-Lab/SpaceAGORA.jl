from __future__ import annotations

from copy import deepcopy
from typing import Any

from orchestrator.app.config import Settings
from orchestrator.app.models import (
    FilePaths,
    MissionIntent,
    OutputOptions,
    SimulationSpec,
    SourceCitation,
)
from orchestrator.app.services.mission_retriever import RetrievedMission


def _deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def _default_payload(intent: MissionIntent, run_name: str) -> dict[str, Any]:
    planet = intent.planet or "Mars"
    return {
        "run_name": run_name,
        "planet": planet,
        "mission_name": intent.mission_name,
        "mission_time_sec": 21600.0,
        "entry_interface_km": 125.0 if planet == "Mars" else 120.0,
        "initial_conditions": {
            "ra_m": 12000e3 if planet == "Mars" else 7000e3,
            "rp_m": 3600e3 if planet == "Mars" else 6800e3,
            "i_deg": 30.0,
            "argument_of_periapsis_deg": 0.0,
            "raan_deg": 0.0,
            "true_anomaly_deg": 180.0,
        },
        "models": {
            "gravity_model": "inverse_squared",
            "density_model": "gram" if planet == "Mars" else "none",
            "thermal_model": "maxwellian",
            "aerodynamic_model": "mach_dependent",
            "wind": planet == "Mars",
            "srp": False,
            "n_bodies": ["Sun"] if planet == "Mars" else [],
        },
        "keplerian": planet == "Earth",
        "orientation_sim": False,
        "num_steps_to_save": 1000,
    }


class SimulationSpecCompiler:
    def __init__(self, settings: Settings) -> None:
        self._settings = settings

    def compile(
        self,
        job_id: str,
        intent: MissionIntent,
        retrieved: RetrievedMission,
        mission_overrides: dict[str, Any],
        output_preferences: dict[str, Any],
    ) -> tuple[SimulationSpec, list[SourceCitation], list[str]]:
        run_name = f"run_{job_id}"
        payload = deepcopy(retrieved.payload) if retrieved.payload else _default_payload(intent, run_name)

        if "run_name" not in payload:
            payload["run_name"] = run_name

        # Carry parsed intent constraints.
        if intent.requested_constraints:
            payload = _deep_update(payload, intent.requested_constraints)

        if mission_overrides:
            payload = _deep_update(payload, mission_overrides)

        # Ensure file paths always come from infrastructure settings.
        payload["file_paths"] = {
            "spice": str(self._settings.spice_path),
            "gram_data": str(self._settings.gram_data_path),
            "gram_py": str(self._settings.gram_py_path),
        }

        # Explicit outputs with user preferences applied.
        output_payload = {
            "results": True,
            "save_csv": True,
            "generate_plots": False,
            "include_log_bundle": True,
        }
        if output_preferences:
            output_payload = _deep_update(output_payload, output_preferences)
        payload["outputs"] = output_payload

        # Build strict typed spec.
        spec = SimulationSpec.model_validate(payload)

        # Normalize mission name from intent when missing.
        if spec.mission_name is None and intent.mission_name is not None:
            spec = spec.model_copy(update={"mission_name": intent.mission_name})

        notes = list(retrieved.retrieval_notes)
        notes.append("compiled_to_simulation_spec")
        return spec, list(retrieved.citations), notes
