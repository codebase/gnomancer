"""Where gnomancer keeps its files.

- Package data (marker panel, haplogroup context, rCRS) ships inside the wheel: `DATA`.
- Downloaded tools (Haplogrep3, PharmCAT, ClinVar) and the optional Python 3.11 side venv live in
  `$GNOMANCER_HOME` (default `~/.gnomancer`). A repo checkout's `tools/` and `.venv-ancestry` are still honoured
  so a development install keeps working.
"""
from __future__ import annotations

import os
from pathlib import Path

PKG = Path(__file__).resolve().parent
DATA = PKG / "data"
REPO = PKG.parent
HOME = Path(os.environ.get("GNOMANCER_HOME") or Path.home() / ".gnomancer")


def tools_dir() -> Path:
    """Directory holding downloaded tools; the first existing candidate, else the home default (created on demand)."""
    for c in (HOME / "tools", REPO / "tools"):
        if c.exists():
            return c
    return HOME / "tools"


def tool(*parts: str) -> Path:
    """Path to a tool file, searching home then repo. Returns the home path if it exists nowhere."""
    for base in (HOME / "tools", REPO / "tools"):
        p = base.joinpath(*parts)
        if p.exists():
            return p
    return (HOME / "tools").joinpath(*parts)


def ancestry_python() -> Path | None:
    for c in (HOME / "venv-ancestry" / "bin" / "python", REPO / ".venv-ancestry" / "bin" / "python"):
        if c.exists():
            return c
    return None


def ancestry_bin(name: str) -> Path | None:
    py = ancestry_python()
    if py and (py.parent / name).exists():
        return py.parent / name
    return None
