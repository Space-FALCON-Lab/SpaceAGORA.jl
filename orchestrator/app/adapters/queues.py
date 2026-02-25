from __future__ import annotations

import json
import queue
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class RunJob:
    job_id: str


class JobQueue(ABC):
    @abstractmethod
    def enqueue(self, job: RunJob) -> None:
        raise NotImplementedError

    @abstractmethod
    def dequeue(self, timeout_sec: float = 1.0) -> RunJob | None:
        raise NotImplementedError


class InMemoryJobQueue(JobQueue):
    def __init__(self) -> None:
        self._queue: queue.Queue[RunJob] = queue.Queue()

    def enqueue(self, job: RunJob) -> None:
        self._queue.put(job)

    def dequeue(self, timeout_sec: float = 1.0) -> RunJob | None:
        try:
            return self._queue.get(timeout=timeout_sec)
        except queue.Empty:
            return None


class SqsJobQueue(JobQueue):
    """AWS SQS-backed queue implementation.

    This adapter is optional and only used when boto3 and AWS credentials are configured.
    """

    def __init__(self, sqs_client: Any, queue_url: str) -> None:
        self._sqs = sqs_client
        self._queue_url = queue_url

    def enqueue(self, job: RunJob) -> None:
        self._sqs.send_message(
            QueueUrl=self._queue_url,
            MessageBody=json.dumps({"job_id": job.job_id}),
        )

    def dequeue(self, timeout_sec: float = 1.0) -> RunJob | None:
        wait_time = int(max(0.0, min(20.0, timeout_sec)))
        response = self._sqs.receive_message(
            QueueUrl=self._queue_url,
            MaxNumberOfMessages=1,
            WaitTimeSeconds=wait_time,
        )
        messages = response.get("Messages", [])
        if not messages:
            return None

        message = messages[0]
        body = json.loads(message["Body"])
        receipt_handle = message["ReceiptHandle"]
        self._sqs.delete_message(QueueUrl=self._queue_url, ReceiptHandle=receipt_handle)
        return RunJob(job_id=body["job_id"])
