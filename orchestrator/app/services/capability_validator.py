from __future__ import annotations

from pathlib import Path

from orchestrator.app.config import Settings
from orchestrator.app.models import CapabilityReport, RunMode, SimulationSpec


class CapabilityValidator:
    def __init__(self, settings: Settings) -> None:
        self._settings = settings

    def validate(self, spec: SimulationSpec, requested_mode: RunMode) -> tuple[SimulationSpec, CapabilityReport]:
        unsupported_reasons: list[str] = []
        notes: list[str] = []
        forced_overrides: dict[str, object] = {}
        reroute_to_async = False

        nearest_spec = spec

        # Rule 1: Earth/Mars only.
        if spec.planet not in {"Earth", "Mars"}:
            unsupported_reasons.append(
                f"planet '{spec.planet}' is not supported in v1; only Earth and Mars are allowed"
            )
            nearest_spec = spec.model_copy(update={"planet": "Mars"})

        # Rule 2: reject near-circular orbits to avoid callback/OE instability.
        ra = spec.initial_conditions.ra_m
        rp = spec.initial_conditions.rp_m
        ecc = abs(ra - rp) / max(ra + rp, 1.0)
        if ecc < 1e-4:
            unsupported_reasons.append(
                "near-circular initial orbit is currently blocked because OE callbacks can fail"
            )

        # Rule 3: force results/save_csv true due current runner behavior.
        if not spec.outputs.results:
            forced_overrides["outputs.results"] = True
            spec = spec.model_copy(update={"outputs": spec.outputs.model_copy(update={"results": True})})
            notes.append("forced outputs.results=true")

        if not spec.outputs.save_csv:
            forced_overrides["outputs.save_csv"] = True
            spec = spec.model_copy(update={"outputs": spec.outputs.model_copy(update={"save_csv": True})})
            notes.append("forced outputs.save_csv=true")

        # Rule 4: required file paths and mission fields.
        self._validate_required_path(Path(spec.file_paths.spice), "SPICE path", unsupported_reasons)
        if spec.models.density_model == "gram":
            self._validate_required_path(Path(spec.file_paths.gram_data), "GRAM data path", unsupported_reasons)
            self._validate_required_path(Path(spec.file_paths.gram_py), "GRAM root path", unsupported_reasons)

        required_spice_files = [
            "pck/pck00011.tpc",
            "lsk/naif0012.tls",
            "spk/planets/de440s.bsp",
        ]
        for relative in required_spice_files:
            full_path = Path(spec.file_paths.spice) / relative
            if not full_path.exists():
                unsupported_reasons.append(f"required SPICE kernel missing: {full_path}")

        if spec.mission_time_sec <= 0:
            unsupported_reasons.append("mission_time_sec must be > 0")

        # Rule 5: reject unsupported models.
        supported_gravity = {"inverse_squared", "inverse_squared_j2", "constant"}
        supported_density = {"gram", "none", "exponential"}
        supported_thermal = {"maxwellian"}
        supported_aero = {"mach_dependent", "constant"}

        if spec.models.gravity_model not in supported_gravity:
            unsupported_reasons.append(
                f"gravity_model '{spec.models.gravity_model}' is unsupported; choose one of {sorted(supported_gravity)}"
            )
        if spec.models.density_model not in supported_density:
            unsupported_reasons.append(
                f"density_model '{spec.models.density_model}' is unsupported; choose one of {sorted(supported_density)}"
            )
        if spec.models.thermal_model not in supported_thermal:
            unsupported_reasons.append(
                f"thermal_model '{spec.models.thermal_model}' is unsupported; choose one of {sorted(supported_thermal)}"
            )
        if spec.models.aerodynamic_model not in supported_aero:
            unsupported_reasons.append(
                f"aerodynamic_model '{spec.models.aerodynamic_model}' is unsupported; choose one of {sorted(supported_aero)}"
            )

        # Rule 6: enforce runtime/resource limits and sync fallback.
        if spec.mission_time_sec > self._settings.max_mission_time_sec:
            unsupported_reasons.append(
                f"mission_time_sec exceeds max allowed ({self._settings.max_mission_time_sec})"
            )

        if requested_mode == RunMode.SYNC and spec.mission_time_sec > self._settings.sync_max_mission_time_sec:
            reroute_to_async = True
            notes.append("sync request rerouted to async due mission_time limit")

        supported = len(unsupported_reasons) == 0
        report = CapabilityReport(
            supported=supported,
            unsupported_reasons=unsupported_reasons,
            nearest_supported_spec=nearest_spec if not supported else None,
            forced_overrides=forced_overrides,
            notes=notes,
            reroute_to_async=reroute_to_async,
        )
        return spec, report

    @staticmethod
    def _validate_required_path(path: Path, label: str, reasons: list[str]) -> None:
        if not path.exists():
            reasons.append(f"{label} does not exist: {path}")
