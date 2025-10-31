"""Centralized logging configuration for OVM-PK."""
import logging
from pathlib import Path
import sys

FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

class LogFilter(logging.Filter):
    """Filter out overly verbose OpenMM logs."""
    def filter(self, record):
        if "OpenMM" in record.name and record.levelno < logging.WARNING:
            return False
        return True

def configure_logging(log_dir: Path = Path("runlogs")):
    """Set up logging with file and console handlers."""
    log_dir.mkdir(exist_ok=True)
    
    # Main handler - writes to file
    main_handler = logging.FileHandler(
        filename=log_dir / "ovmpk.log",
        mode="a"
    )
    main_handler.setFormatter(logging.Formatter(FORMAT, DATE_FORMAT))
    main_handler.addFilter(LogFilter())
    
    # Console handler - only warnings+
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(logging.Formatter(FORMAT, DATE_FORMAT))
    
    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        handlers=[main_handler, console_handler],
        force=True
    )
    
    # Special cases
    logging.getLogger("openmm").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
