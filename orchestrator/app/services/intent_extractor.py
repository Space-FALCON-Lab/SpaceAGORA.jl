from __future__ import annotations

import re
from abc import ABC, abstractmethod

from orchestrator.app.models import MissionIntent


class IntentExtractor(ABC):
    @abstractmethod
    def extract(self, prompt: str) -> MissionIntent:
        raise NotImplementedError


class RuleBasedIntentExtractor(IntentExtractor):
    """Deterministic intent extraction for pilot phase.

    In production, this class can be replaced by an LLM-backed extractor while
    preserving the same MissionIntent contract.
    """

    _MARS_ODYSSEY_PATTERNS = (
        r"mars\s+odyssey",
        r"odyssey\s+mission",
        r"replicat(e|es|ing)\s+.*odyssey",
    )

    def extract(self, prompt: str) -> MissionIntent:
        text = prompt.strip()
        lower = text.lower()

        planet = None
        if "mars" in lower:
            planet = "Mars"
        elif "earth" in lower:
            planet = "Earth"

        mission_name = None
        replicate = False
        for pattern in self._MARS_ODYSSEY_PATTERNS:
            if re.search(pattern, lower):
                mission_name = "Mars Odyssey"
                replicate = True
                planet = "Mars"
                break

        constraints: dict[str, object] = {}

        # Basic mission-time parsing from free text.
        match_hours = re.search(r"(\d+(?:\.\d+)?)\s*hours?", lower)
        match_minutes = re.search(r"(\d+(?:\.\d+)?)\s*minutes?", lower)
        match_days = re.search(r"(\d+(?:\.\d+)?)\s*days?", lower)

        if match_days:
            constraints["mission_time_sec"] = float(match_days.group(1)) * 86400.0
        elif match_hours:
            constraints["mission_time_sec"] = float(match_hours.group(1)) * 3600.0
        elif match_minutes:
            constraints["mission_time_sec"] = float(match_minutes.group(1)) * 60.0

        objective = "simulate_mission"
        if "initial condition" in lower or "initial conditions" in lower:
            objective = "simulate_with_custom_initial_conditions"

        return MissionIntent(
            objective=objective,
            mission_name=mission_name,
            planet=planet,
            replicate_mission=replicate,
            requested_constraints=constraints,
        )
