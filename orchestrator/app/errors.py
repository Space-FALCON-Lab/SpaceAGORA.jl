class OrchestratorError(Exception):
    """Base error for orchestrator failures."""


class NotFoundError(OrchestratorError):
    """Raised when a requested entity does not exist."""


class UnsupportedRequestError(OrchestratorError):
    """Raised when a request cannot be supported by current capabilities."""


class ExecutionError(OrchestratorError):
    """Raised when simulation execution fails."""
