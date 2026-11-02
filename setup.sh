#!/usr/bin/env bash
# Developer setup from a checkout: venv + editable install, then `gnomancer setup` for tools and references.
# End users: pipx install gnomancer && gnomancer setup   (or the installer at https://gnomancer.com/install.sh)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
[ -d .venv ] || $PY -m venv .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -e ".[mcp,dev]"
echo "python deps ok"
.venv/bin/gnomancer setup "$@"
echo
echo "Run:  .venv/bin/gnomancer /path/to/raw-file.txt [--ancestry]"
