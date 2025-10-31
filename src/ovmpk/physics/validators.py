"""Calibration guardrails for physics modules."""

class CalibrationForbiddenError(RuntimeError):
    """Raised when empirical adjustments are detected."""
    def __init__(self, param_name: str):
        super().__init__(
            f"Calibration parameter '{param_name}' is forbidden. "
            "Use physical models (MCPB/QMMM) instead."
        )

def validate_no_calibration(params: dict):
    """Check for empirical correction terms."""
    forbidden = {'bonus', 'offset', 'scale_factor', 'calibration'}
    for param in params:
        if any(f in param.lower() for f in forbidden):
            raise CalibrationForbiddenError(param)
