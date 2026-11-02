#!/usr/bin/env bash
# gnomancer installer for macOS and Linux.   Usage:  curl -fsSL https://gnomancer.com/install.sh | bash
# Installs the PyPI package into ~/.gnomancer/venv and downloads its tools. Re-running upgrades it.
# Nothing is uploaded; see https://gnomancer.com/#privacy
set -euo pipefail
HOME_DIR="${GNOMANCER_HOME:-$HOME/.gnomancer}"
say() { printf '\033[1m%s\033[0m\n' "$*"; }
need() { command -v "$1" >/dev/null 2>&1; }

say "gnomancer installer"
# A usable Python: 3.10+, and able to load its XML module (Homebrew's Python bottles on macOS 26.2 cannot, which breaks pip).
healthy() { "$1" -c 'import sys, platform, pyexpat; sys.exit(0 if sys.version_info >= (3,10) and platform.mac_ver()[0] != "" or sys.platform != "darwin" else 1)' 2>/dev/null; }
PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3; do
  if need "$c" && healthy "$(command -v "$c")"; then PY="$(command -v "$c")"; break; fi
done
if [ -z "$PY" ] && need uv; then
  say "No healthy system Python found; using uv to fetch a standalone CPython 3.13"
  uv python install 3.13 -q && PY="$(uv python find 3.13)"
fi
if [ -z "$PY" ]; then
  echo "No working Python 3.10+ found (a Python that is installed but cannot import its XML module counts as not working)."
  echo "Easiest fix on any platform: install uv (https://docs.astral.sh/uv/) and re-run this script; it will fetch a standalone Python."
  if [ "$(uname)" = "Darwin" ]; then echo "macOS: install Homebrew (https://brew.sh), then:  brew install python@3.12 openjdk";
  else echo "Debian/Ubuntu:  sudo apt install python3 python3-venv default-jre     Fedora:  sudo dnf install python3 java-21-openjdk"; fi
  exit 1
fi
if ! need java && [ ! -x /opt/homebrew/opt/openjdk/bin/java ] && [ ! -x /usr/local/opt/openjdk/bin/java ]; then
  echo "Note: Java 11+ not found. Mitochondrial haplogroups and PharmCAT will be skipped until it is installed"
  echo "      (macOS: brew install openjdk   Debian/Ubuntu: sudo apt install default-jre). Everything else works."
fi

mkdir -p "$HOME_DIR"
if [ ! -x "$HOME_DIR/venv/bin/python" ]; then say "Creating $HOME_DIR/venv with $PY"; "$PY" -m venv "$HOME_DIR/venv"; fi
say "Installing gnomancer from PyPI"
"$HOME_DIR/venv/bin/pip" install -q --upgrade pip
"$HOME_DIR/venv/bin/pip" install -q --upgrade "gnomancer[mcp]"
say "Downloading tools and reference files (once; ~250 MB, mostly ClinVar)"
GNOMANCER_HOME="$HOME_DIR" "$HOME_DIR/venv/bin/gnomancer" setup
mkdir -p "$HOME/.local/bin"
ln -sf "$HOME_DIR/venv/bin/gnomancer" "$HOME/.local/bin/gnomancer"
say "Done."
echo
echo "Run it on your raw data file (how to download that: https://gnomancer.com/get-your-data/):"
echo
echo "    $HOME_DIR/venv/bin/gnomancer ~/Downloads/AncestryDNA.zip"
echo
case ":$PATH:" in *":$HOME/.local/bin:"*) echo "or simply:   gnomancer ~/Downloads/AncestryDNA.zip";; *) echo "Add ~/.local/bin to your PATH to run it as just 'gnomancer'.";; esac
