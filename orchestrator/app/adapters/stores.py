from __future__ import annotations

import json
import threading
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Callable

from orchestrator.app.errors import NotFoundError
from orchestrator.app.models import RunRecord


class RunStore(ABC):
    @abstractmethod
    def put(self, record: RunRecord) -> None:
        raise NotImplementedError

    @abstractmethod
    def get(self, job_id: str) -> RunRecord:
        raise NotImplementedError

    @abstractmethod
    def update(self, job_id: str, updater: Callable[[RunRecord], RunRecord]) -> RunRecord:
        raise NotImplementedError


class FileRunStore(RunStore):
    def __init__(self, runs_dir: Path) -> None:
        self._runs_dir = runs_dir
        self._runs_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _path_for(self, job_id: str) -> Path:
        return self._runs_dir / f"{job_id}.json"

    def put(self, record: RunRecord) -> None:
        path = self._path_for(record.job_id)
        with self._lock:
            path.write_text(record.model_dump_json(indent=2), encoding="utf-8")

    def get(self, job_id: str) -> RunRecord:
        path = self._path_for(job_id)
        if not path.exists():
            raise NotFoundError(f"Run '{job_id}' not found")
        payload = json.loads(path.read_text(encoding="utf-8"))
        return RunRecord.model_validate(payload)

    def update(self, job_id: str, updater: Callable[[RunRecord], RunRecord]) -> RunRecord:
        with self._lock:
            current = self.get(job_id)
            updated = updater(current)
            self.put(updated)
            return updated


class DynamoRunStore(RunStore):
    """AWS DynamoDB-backed run store.

    This adapter is optional and only used when boto3 and table config are provided.
    """

    def __init__(self, table: Any) -> None:
        self._table = table

    def put(self, record: RunRecord) -> None:
        self._table.put_item(Item=record.model_dump(mode="json"))

    def get(self, job_id: str) -> RunRecord:
        response = self._table.get_item(Key={"job_id": job_id})
        item = response.get("Item")
        if not item:
            raise NotFoundError(f"Run '{job_id}' not found")
        return RunRecord.model_validate(item)

    def update(self, job_id: str, updater: Callable[[RunRecord], RunRecord]) -> RunRecord:
        current = self.get(job_id)
        updated = updater(current)
        self.put(updated)
        return updated
