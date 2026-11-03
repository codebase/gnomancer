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


JAVA_HOMES = ("/opt/homebrew/opt/openjdk", "/usr/local/opt/openjdk", "/usr/lib/jvm/default-java")


def java_bin() -> str | None:
    """Path to a `java` executable: $JAVA_HOME, then the usual Homebrew/Debian locations, then PATH."""
    import shutil
    for jh in (os.environ.get("JAVA_HOME"), *JAVA_HOMES):
        if jh and Path(jh, "bin", "java").exists():
            return str(Path(jh, "bin", "java"))
    return shutil.which("java")


def java_env() -> dict:
    """A copy of the environment with JAVA_HOME/PATH set so a subprocess finds the same `java` as java_bin()."""
    env = dict(os.environ)
    if "JAVA_HOME" not in env:
        for jh in JAVA_HOMES:
            if Path(jh, "bin", "java").exists():
                env["JAVA_HOME"] = jh
                env["PATH"] = f"{jh}/bin:" + env.get("PATH", "")
                break
    return env


def local_bin(name: str) -> str | None:
    """An executable installed next to the running interpreter (same venv), else on PATH."""
    import shutil, sys
    cand = Path(sys.executable).parent / name
    return str(cand) if cand.exists() else shutil.which(name)


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
