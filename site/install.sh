#!/usr/bin/env bash
# gnomancer installer for macOS and Linux.   Usage:  curl -fsSL https://gnomancer.com/install.sh | bash
# Installs the PyPI package into ~/.gnomancer/venv and downloads its tools. Re-running upgrades it.
# Nothing is uploaded; see https://gnomancer.com/#privacy
set -euo pipefail
HOME_DIR="${GNOMANCER_HOME:-$HOME/.gnomancer}"
say() { printf '\033[1m%s\033[0m\n' "$*"; }
need() { command -v "$1" >/dev/null 2>&1; }

say "gnomancer installer"
# A usable Python: 3.10+, able to create a venv, and able to load its XML module
# (Homebrew's Python bottles on macOS 26.2 cannot, which breaks pip).
healthy() { "$1" -c 'import sys, venv, pyexpat; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; }
PY=""
for c in python3.13 python3.12 python3.11 python3.10 python3; do
  if need "$c" && healthy "$(command -v "$c")"; then PY="$(command -v "$c")"; break; fi
done
if [ -z "$PY" ] && need uv; then
  say "No healthy system Python found; using uv to fetch a standalone CPython 3.13"
  uv python install 3.13 -q && PY="$(uv python find 3.13)"
fi
if [ -z "$PY" ]; then
  echo "No working Python 3.10+ found (a Python that cannot create a venv or import its XML module counts as not working)."
  echo "Easiest fix on any platform: install uv (https://docs.astral.sh/uv/) and re-run this script; it fetches a standalone Python."
  if [ "$(uname)" = "Darwin" ]; then echo "Or download Python from https://www.python.org/downloads/ (the python.org build works where Homebrew's currently does not).";
  else echo "Debian/Ubuntu:  sudo apt install python3 python3-venv     Fedora:  sudo dnf install python3"; fi
  exit 1
fi

# Java 17+ is needed only for the mtDNA haplogroup (Haplogrep3) and PharmCAT steps.
JAVA=""
for j in "${JAVA_HOME:-}/bin/java" /opt/homebrew/opt/openjdk/bin/java /usr/local/opt/openjdk/bin/java "$(command -v java 2>/dev/null || true)"; do
  [ -n "$j" ] && [ -x "$j" ] && { JAVA="$j"; break; }
done
if [ -n "$JAVA" ]; then
  MAJOR="$("$JAVA" -version 2>&1 | head -1 | sed -E 's/.*"([0-9]+)(\.[0-9]+)*.*/\1/')"
  case "$MAJOR" in ''|*[!0-9]*) MAJOR=0;; esac
  [ "$MAJOR" -eq 1 ] && MAJOR=8
  if [ "$MAJOR" -lt 17 ]; then echo "Note: Java $MAJOR found; the mtDNA haplogroup and PharmCAT steps need Java 17+ and will be skipped until it is upgraded."; fi
else
  echo "Note: Java not found. The mtDNA haplogroup and PharmCAT steps will be skipped until Java 17+ is installed"
  echo "      (macOS: brew install openjdk   Debian/Ubuntu: sudo apt install openjdk-17-jre   Fedora: sudo dnf install java-21-openjdk). Everything else works."
fi

mkdir -p "$HOME_DIR"
# Repair a venv whose base Python has since broken (common after a Homebrew python upgrade).
if [ -x "$HOME_DIR/venv/bin/python" ] && ! "$HOME_DIR/venv/bin/python" -c 'import pyexpat' 2>/dev/null; then
  say "Existing venv is broken (its Python can no longer import pyexpat); recreating it"
  rm -rf "$HOME_DIR/venv"
fi
if [ ! -x "$HOME_DIR/venv/bin/python" ]; then say "Creating $HOME_DIR/venv with $PY"; "$PY" -m venv "$HOME_DIR/venv"; fi
say "Installing gnomancer from PyPI"
"$HOME_DIR/venv/bin/python" -m pip install -q --upgrade pip
"$HOME_DIR/venv/bin/python" -m pip install -q --upgrade gnomancer
say "Downloading tools and reference files (once; ~250 MB, mostly ClinVar)"
GNOMANCER_HOME="$HOME_DIR" "$HOME_DIR/venv/bin/gnomancer" setup
mkdir -p "$HOME/.local/bin"
ln -sf "$HOME_DIR/venv/bin/gnomancer" "$HOME/.local/bin/gnomancer"
say "Done."
echo
echo "Run it on your raw data file (how to download that: https://gnomancer.com/get-your-data/):"
echo
echo "    gnomancer ~/Downloads/AncestryDNA.zip"
echo
case ":$PATH:" in
  *":$HOME/.local/bin:"*) ;;
  *)
    RC="$HOME/.bashrc"; case "${SHELL:-}" in */zsh) RC="$HOME/.zshrc";; esac
    echo "~/.local/bin is not on your PATH yet. Add it once, then open a new terminal:"
    echo
    echo "    echo 'export PATH=\"\$HOME/.local/bin:\$PATH\"' >> $RC"
    echo
    echo "Until then, the full path works:  $HOME_DIR/venv/bin/gnomancer ~/Downloads/AncestryDNA.zip";;
esac
