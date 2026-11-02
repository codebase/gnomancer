import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="session")
def demo_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("demo")
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_demo.py"), str(d)], check=True, capture_output=True)
    return d
