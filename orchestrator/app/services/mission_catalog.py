from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from orchestrator.app.models import SourceCitation


@dataclass(frozen=True)
class MissionTemplate:
    mission_id: str
    mission_name: str
    payload: dict[str, Any]
    citations: list[SourceCitation]


class CuratedMissionCatalog:
    """Curated mission templates for deterministic phase-1 runs."""

    def __init__(self) -> None:
        self._templates = {
            "mars_odyssey": MissionTemplate(
                mission_id="mars_odyssey",
                mission_name="Mars Odyssey",
                payload={
                    "run_name": "mars_odyssey",
                    "planet": "Mars",
                    "mission_name": "Mars Odyssey",
                    "mission_time_sec": 86400.0,
                    "entry_interface_km": 125.0,
                    "initial_conditions": {
                        "ra_m": 28559.615e3,
                        "rp_m": 3396.0e3 + 87.0e3,
                        "i_deg": 93.522,
                        "argument_of_periapsis_deg": 109.7454,
                        "raan_deg": 28.1517,
                        "true_anomaly_deg": 180.0,
                    },
                    "models": {
                        "gravity_model": "inverse_squared",
                        "density_model": "gram",
                        "thermal_model": "maxwellian",
                        "aerodynamic_model": "mach_dependent",
                        "wind": True,
                        "srp": False,
                        "n_bodies": ["Sun"],
                    },
                    "keplerian": False,
                    "orientation_sim": False,
                    "num_steps_to_save": 1000,
                },
                citations=[
                    SourceCitation(
                        label="SpaceAGORA Mars Odyssey Example",
                        url="https://github.com/Space-FALCON-Lab/ABTS.jl/blob/main/src/examples/AGORA_Odyssey.jl",
                    ),
                    SourceCitation(
                        label="NASA Mars Odyssey Mission Overview",
                        url="https://science.nasa.gov/mission/odyssey/",
                    ),
                    SourceCitation(
                        label="JPL Mars Odyssey Mission Profile",
                        url="https://www.jpl.nasa.gov/missions/2001-mars-odyssey",
                    ),
                ],
            ),
            "generic_mars": MissionTemplate(
                mission_id="generic_mars",
                mission_name="Generic Mars Orbiter",
                payload={
                    "run_name": "generic_mars",
                    "planet": "Mars",
                    "mission_time_sec": 21600.0,
                    "entry_interface_km": 125.0,
                    "initial_conditions": {
                        "ra_m": 12000e3,
                        "rp_m": 3600e3,
                        "i_deg": 30.0,
                        "argument_of_periapsis_deg": 0.0,
                        "raan_deg": 0.0,
                        "true_anomaly_deg": 180.0,
                    },
                    "models": {
                        "gravity_model": "inverse_squared",
                        "density_model": "gram",
                        "thermal_model": "maxwellian",
                        "aerodynamic_model": "mach_dependent",
                        "wind": True,
                        "srp": False,
                        "n_bodies": ["Sun"],
                    },
                    "keplerian": False,
                    "orientation_sim": False,
                    "num_steps_to_save": 1000,
                },
                citations=[
                    SourceCitation(
                        label="SpaceAGORA Documentation",
                        url="https://github.com/Space-FALCON-Lab/ABTS.jl",
                    )
                ],
            ),
            "generic_earth": MissionTemplate(
                mission_id="generic_earth",
                mission_name="Generic Earth Orbiter",
                payload={
                    "run_name": "generic_earth",
                    "planet": "Earth",
                    "mission_time_sec": 21600.0,
                    "entry_interface_km": 120.0,
                    "initial_conditions": {
                        "ra_m": 7000e3,
                        "rp_m": 6800e3,
                        "i_deg": 28.5,
                        "argument_of_periapsis_deg": 0.0,
                        "raan_deg": 0.0,
                        "true_anomaly_deg": 180.0,
                    },
                    "models": {
                        "gravity_model": "inverse_squared",
                        "density_model": "none",
                        "thermal_model": "maxwellian",
                        "aerodynamic_model": "mach_dependent",
                        "wind": False,
                        "srp": False,
                        "n_bodies": [],
                    },
                    "keplerian": True,
                    "orientation_sim": False,
                    "num_steps_to_save": 1000,
                },
                citations=[
                    SourceCitation(
                        label="SpaceAGORA Earth Example",
                        url="https://github.com/Space-FALCON-Lab/ABTS.jl/blob/main/src/examples/AGORA_Earth.jl",
                    )
                ],
            ),
        }

    def get(self, mission_id: str) -> MissionTemplate | None:
        return self._templates.get(mission_id)

    def resolve_id(self, mission_name: str | None, planet: str | None) -> str | None:
        if mission_name is None:
            if planet == "Mars":
                return "generic_mars"
            if planet == "Earth":
                return "generic_earth"
            return None

        key = mission_name.strip().lower()
        if key in {"mars odyssey", "odyssey", "2001 mars odyssey"}:
            return "mars_odyssey"
        if planet == "Mars":
            return "generic_mars"
        if planet == "Earth":
            return "generic_earth"
        return None
