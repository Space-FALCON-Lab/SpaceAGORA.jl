from __future__ import annotations

import mimetypes
import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from orchestrator.app.models import ArtifactRef


class ArtifactStore(ABC):
    @abstractmethod
    def save_file(self, job_id: str, source: Path, artifact_id: str, kind: str) -> ArtifactRef:
        raise NotImplementedError

    @abstractmethod
    def resolve_path(self, job_id: str, artifact_id: str) -> Path:
        raise NotImplementedError


class LocalArtifactStore(ArtifactStore):
    def __init__(self, root: Path) -> None:
        self._root = root
        self._root.mkdir(parents=True, exist_ok=True)

    def _job_dir(self, job_id: str) -> Path:
        job_dir = self._root / job_id
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    def save_file(self, job_id: str, source: Path, artifact_id: str, kind: str) -> ArtifactRef:
        ext = source.suffix
        target = self._job_dir(job_id) / f"{artifact_id}{ext}"
        shutil.copy2(source, target)
        content_type = mimetypes.guess_type(str(target))[0] or "application/octet-stream"
        return ArtifactRef(
            artifact_id=artifact_id,
            kind=kind,
            path=str(target),
            content_type=content_type,
            size_bytes=target.stat().st_size,
        )

    def resolve_path(self, job_id: str, artifact_id: str) -> Path:
        job_dir = self._job_dir(job_id)
        candidates = list(job_dir.glob(f"{artifact_id}.*"))
        if not candidates:
            raise FileNotFoundError(f"Artifact '{artifact_id}' for run '{job_id}' not found")
        return candidates[0]


class S3ArtifactStore(ArtifactStore):
    """AWS S3-backed artifact store.

    For direct downloads, API can return pre-signed URLs instead of local paths.
    """

    def __init__(self, s3_client: Any, bucket: str, key_prefix: str = "runs") -> None:
        self._s3 = s3_client
        self._bucket = bucket
        self._prefix = key_prefix.strip("/")

    def save_file(self, job_id: str, source: Path, artifact_id: str, kind: str) -> ArtifactRef:
        key = f"{self._prefix}/{job_id}/{artifact_id}{source.suffix}"
        content_type = mimetypes.guess_type(str(source))[0] or "application/octet-stream"
        self._s3.upload_file(
            Filename=str(source),
            Bucket=self._bucket,
            Key=key,
            ExtraArgs={"ContentType": content_type},
        )
        return ArtifactRef(
            artifact_id=artifact_id,
            kind=kind,
            path=f"s3://{self._bucket}/{key}",
            content_type=content_type,
            size_bytes=source.stat().st_size,
        )

    def resolve_path(self, job_id: str, artifact_id: str) -> Path:
        raise NotImplementedError("S3ArtifactStore does not provide local paths")

    def presigned_url(self, artifact: ArtifactRef, expires_sec: int = 3600) -> str:
        path = artifact.path.replace("s3://", "", 1)
        bucket, key = path.split("/", 1)
        return self._s3.generate_presigned_url(
            ClientMethod="get_object",
            Params={"Bucket": bucket, "Key": key},
            ExpiresIn=expires_sec,
        )
