from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from orchestrator.app.models import ArtifactRef, ResultSummary, SimulationSpec, SourceCitation


def summarize_result(
    spec: SimulationSpec,
    result_csv: Path | None,
    citations: list[SourceCitation],
    artifacts: list[ArtifactRef],
    extra_caveats: list[str] | None = None,
) -> ResultSummary:
    caveats = list(extra_caveats or [])
    metrics: dict[str, Any] = {
        "planet": spec.planet,
        "mission_name": spec.mission_name,
        "mission_time_sec": spec.mission_time_sec,
    }

    if result_csv is None or not result_csv.exists():
        caveats.append("simulation_results.csv was not produced by the runner")
        summary_text = (
            f"Run completed without a parsed trajectory table for {spec.planet}. "
            "Check logs/artifacts for runtime details."
        )
        return ResultSummary(
            summary_text=summary_text,
            key_metrics=metrics,
            caveats=caveats,
            citations=citations,
            artifacts=artifacts,
        )

    rows = []
    with result_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)

    metrics["rows"] = len(rows)
    if rows:
        first = rows[0]
        last = rows[-1]
        metrics["start_time_s"] = float(first.get("time", 0.0))
        metrics["end_time_s"] = float(last.get("time", 0.0))

        try:
            x = float(last.get("sc1_pos_1", 0.0))
            y = float(last.get("sc1_pos_2", 0.0))
            z = float(last.get("sc1_pos_3", 0.0))
            vx = float(last.get("sc1_vel_1", 0.0))
            vy = float(last.get("sc1_vel_2", 0.0))
            vz = float(last.get("sc1_vel_3", 0.0))
            metrics["final_radius_m"] = math.sqrt(x * x + y * y + z * z)
            metrics["final_speed_mps"] = math.sqrt(vx * vx + vy * vy + vz * vz)
        except ValueError:
            caveats.append("unable to parse final state numerically from simulation_results.csv")

    summary_text = (
        f"Simulation finished for {spec.planet}"
        + (f" ({spec.mission_name})" if spec.mission_name else "")
        + f" with {metrics.get('rows', 0)} saved state rows."
    )

    if not citations:
        caveats.append("no mission data citations were attached")

    return ResultSummary(
        summary_text=summary_text,
        key_metrics=metrics,
        caveats=caveats,
        citations=citations,
        artifacts=artifacts,
    )
