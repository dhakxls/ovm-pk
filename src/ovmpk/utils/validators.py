"""Validation utilities for physics modules."""

class CalibrationForbiddenError(RuntimeError):
    """Raised when empirical calibration is attempted."""
    def __init__(self, param_name: str):
        super().__init__(f"Empirical parameter '{param_name}' is forbidden. Use physical models instead")

def validate_no_calibration(params: dict):
    """Check for disallowed calibration parameters."""
    forbidden_terms = {'bonus', 'offset', 'scale_factor', 'calibration'}
    for param in params:
        if any(term in param.lower() for term in forbidden_terms):
            raise CalibrationForbiddenError(param)
