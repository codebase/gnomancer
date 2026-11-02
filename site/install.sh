#!/usr/bin/env bash
# gnomancer installer for macOS and Linux.   Usage:  curl -fsSL https://gnomancer.com/install.sh | bash
# Installs into ~/gnomancer. Re-running updates it. Nothing is uploaded; see https://gnomancer.com/#privacy
set -euo pipefail
DIR="${GNOMANCER_DIR:-$HOME/gnomancer}"
say() { printf '\033[1m%s\033[0m\n' "$*"; }
need() { command -v "$1" >/dev/null 2>&1; }

say "gnomancer installer"
missing=()
need git || missing+=(git)
if need python3 && python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3,10) else 1)'; then :; else missing+=("python3 (3.10 or newer)"); fi
if [ ${#missing[@]} -gt 0 ]; then
  echo "Missing: ${missing[*]}"
  if [ "$(uname)" = "Darwin" ]; then
    echo "On macOS the easiest fix is Homebrew (https://brew.sh), then:  brew install git python@3.12 openjdk"
  else
    echo "On Debian/Ubuntu:  sudo apt install git python3 python3-venv default-jre     On Fedora:  sudo dnf install git python3 java-21-openjdk"
  fi
  exit 1
fi
if ! need java && [ ! -x /opt/homebrew/opt/openjdk/bin/java ] && [ ! -x /usr/local/opt/openjdk/bin/java ]; then
  echo "Note: Java 11+ was not found. The haplogroup (Haplogrep3) and PharmCAT steps will be skipped until it is installed"
  echo "      (macOS: brew install openjdk   Debian/Ubuntu: sudo apt install default-jre). Everything else works without it."
fi

if [ -d "$DIR/.git" ]; then
  say "Updating $DIR"; git -C "$DIR" pull -q --ff-only
else
  say "Downloading to $DIR"; git clone -q https://github.com/codebase/gnomancer "$DIR"
fi
cd "$DIR"
say "Setting up (this downloads the open-source tools and a ~200 MB ClinVar file once)"
./setup.sh
mkdir -p "$HOME/.local/bin"
ln -sf "$DIR/.venv/bin/gnomancer" "$HOME/.local/bin/gnomancer"
say "Done."
echo
echo "Run it on your raw data file (see https://gnomancer.com/get-your-data/ for how to download that):"
echo
echo "    $DIR/.venv/bin/gnomancer ~/Downloads/AncestryDNA.zip"
echo
case ":$PATH:" in *":$HOME/.local/bin:"*) echo "or simply:   gnomancer ~/Downloads/AncestryDNA.zip";; *) echo "Add ~/.local/bin to your PATH to run it as just 'gnomancer'.";; esac
