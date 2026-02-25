from __future__ import annotations

import threading
import time
from typing import Optional

from orchestrator.app.services.run_orchestrator import RunOrchestrator


class BackgroundWorker:
    def __init__(self, orchestrator: RunOrchestrator, poll_interval_sec: float = 0.5) -> None:
        self._orchestrator = orchestrator
        self._poll_interval_sec = poll_interval_sec
        self._thread: Optional[threading.Thread] = None
        self._stop_event = threading.Event()

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._run_loop, name="orchestrator-worker", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        if self._thread is not None:
            self._thread.join(timeout=3.0)

    def _run_loop(self) -> None:
        while not self._stop_event.is_set():
            try:
                record = self._orchestrator.execute_next_queued_job()
                if record is None:
                    time.sleep(self._poll_interval_sec)
            except Exception:
                # Keep worker alive for subsequent jobs.
                time.sleep(self._poll_interval_sec)
