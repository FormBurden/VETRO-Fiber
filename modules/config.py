# modules/config.py
# Central place to define the project data directory.
# Per project rules, ALWAYS reference this as: modules.config.DATA_DIR

import os
from pathlib import Path

def _default_data_dir() -> Path:
    """
    Default to a 'data' folder at the repo root (one level above modules/).
    """
    modules_dir = Path(__file__).resolve().parent
    repo_root = modules_dir.parent
    return repo_root / "data"

# Allow override via environment variable (optional).
# Example:
#   Windows PowerShell:  $env:VETRO_DATA_DIR = "D:\VETRO\data"
#   Linux/macOS bash:    export VETRO_DATA_DIR="$HOME/vetro_data"
_env = os.environ.get("VETRO_DATA_DIR")

if _env:
    _dir = Path(os.path.expanduser(_env)).resolve()
else:
    _dir = _default_data_dir().resolve()

# Ensure the directory exists.
_dir.mkdir(parents=True, exist_ok=True)

# Export as a string path for easy use with os.path.join and Path alike.
DATA_DIR = str(_dir)

__all__ = ["DATA_DIR"]
