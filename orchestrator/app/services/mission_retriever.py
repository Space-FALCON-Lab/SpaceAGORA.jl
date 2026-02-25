from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Optional

from orchestrator.app.models import MissionIntent, SourceCitation
from orchestrator.app.services.mission_catalog import CuratedMissionCatalog


@dataclass
class RetrievedMission:
    mission_id: Optional[str] = None
    payload: dict[str, Any] = field(default_factory=dict)
    citations: list[SourceCitation] = field(default_factory=list)
    retrieval_notes: list[str] = field(default_factory=list)


class MissionRetriever:
    """Hybrid mission retriever: curated catalog first, optional web fallback."""

    def __init__(self, catalog: CuratedMissionCatalog, enable_web_fallback: bool) -> None:
        self._catalog = catalog
        self._enable_web_fallback = enable_web_fallback

    def retrieve(self, intent: MissionIntent) -> RetrievedMission:
        mission_id = self._catalog.resolve_id(intent.mission_name, intent.planet)
        if mission_id:
            template = self._catalog.get(mission_id)
            if template is not None:
                return RetrievedMission(
                    mission_id=mission_id,
                    payload=dict(template.payload),
                    citations=list(template.citations),
                    retrieval_notes=["resolved_from_curated_catalog"],
                )

        result = RetrievedMission()
        result.retrieval_notes.append("curated_catalog_miss")
        if not self._enable_web_fallback:
            result.retrieval_notes.append("web_fallback_disabled")
            return result

        if intent.mission_name:
            web_payload, web_citations = self._retrieve_from_web(intent)
            result.payload.update(web_payload)
            result.citations.extend(web_citations)
            if web_payload:
                result.retrieval_notes.append("filled_from_web_retrieval")
            else:
                result.retrieval_notes.append("web_retrieval_no_structured_fields")

        return result

    def _retrieve_from_web(self, intent: MissionIntent) -> tuple[dict[str, Any], list[SourceCitation]]:
        """Best-effort retrieval stub for pilot phase.

        This intentionally returns conservative structured fields and citations,
        and avoids fragile scraping logic.
        """

        mission_name = intent.mission_name or "Unknown Mission"
        query = urllib.parse.quote_plus(mission_name)
        nasa_search = f"https://api.nasa.gov/planetary/apod?api_key=DEMO_KEY"  # Reachability probe
        jpl_search = f"https://www.jpl.nasa.gov/search?q={query}"

        citations: list[SourceCitation] = [
            SourceCitation(
                label="NASA Search",
                url=f"https://www.nasa.gov/search/?q={query}",
                retrieved_at=datetime.now(timezone.utc),
            ),
            SourceCitation(
                label="JPL Search",
                url=jpl_search,
                retrieved_at=datetime.now(timezone.utc),
            ),
        ]

        # Optional lightweight probe to confirm outbound access for web retrieval.
        try:
            req = urllib.request.Request(nasa_search, method="GET")
            with urllib.request.urlopen(req, timeout=3) as response:
                _ = response.read(32)
        except (urllib.error.URLError, TimeoutError, ValueError):
            pass

        payload: dict[str, Any] = {
            "mission_name": mission_name,
        }
        if intent.planet:
            payload["planet"] = intent.planet

        return payload, citations
