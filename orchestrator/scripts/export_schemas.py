from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from orchestrator.app.models import (
    CapabilityReport,
    MissionIntent,
    ResultSummary,
    RunRecord,
    SimulationSpec,
)


def main() -> None:
    out_dir = Path(__file__).resolve().parents[1] / "schemas"
    out_dir.mkdir(parents=True, exist_ok=True)

    schemas = {
        "MissionIntent": MissionIntent.model_json_schema(),
        "SimulationSpec": SimulationSpec.model_json_schema(),
        "CapabilityReport": CapabilityReport.model_json_schema(),
        "RunRecord": RunRecord.model_json_schema(),
        "ResultSummary": ResultSummary.model_json_schema(),
    }

    for name, schema in schemas.items():
        path = out_dir / f"{name}.schema.json"
        path.write_text(json.dumps(schema, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
