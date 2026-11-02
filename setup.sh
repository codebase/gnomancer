#!/usr/bin/env bash
# One-time setup: Python venv with snps + yhaplo, and Haplogrep3 (needs Java 11+).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
[ -d .venv ] || $PY -m venv .venv
.venv/bin/pip install -q --upgrade pip
.venv/bin/pip install -q -e ".[mcp]" -r requirements.txt
echo "python deps ok"

# Haplogrep3 (mtDNA). The "linux" bundle is a jar + launcher and runs on macOS too.
HG_VER=${HAPLOGREP3_VERSION:-3.3.2}
if [ ! -x tools/haplogrep3/haplogrep3 ]; then
  mkdir -p tools/haplogrep3
  curl -sL "https://github.com/genepi/haplogrep3/releases/download/v${HG_VER}/haplogrep3-${HG_VER}-linux.tar.gz" \
    | tar xz -C tools/haplogrep3
fi
# Java: prefer Homebrew OpenJDK on macOS if the system java is too old.
for jh in /opt/homebrew/opt/openjdk /usr/local/opt/openjdk /usr/lib/jvm/default-java; do
  if [ -x "$jh/bin/java" ]; then export JAVA_HOME="$jh"; export PATH="$jh/bin:$PATH"; break; fi
done
if java -version 2>&1 | head -1 | grep -qE '"(1\.[0-9]|9|10)\.'; then
  echo "WARNING: java is older than 11; Haplogrep3 will not run. Install a JDK (e.g. brew install openjdk)." >&2
else
  echo "haplogrep3 ok ($(java -version 2>&1 | head -1))"
fi

# ClinVar (GRCh37 VCF, ~200 MB) for the whole-chip pathogenic-variant scan.
if [ ! -s tools/clinvar/clinvar_GRCh37.vcf.gz ]; then
  mkdir -p tools/clinvar
  echo "downloading ClinVar GRCh37 (~200 MB) ..."
  curl -sL -o tools/clinvar/clinvar_GRCh37.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
fi
echo "clinvar ok ($(ls -la tools/clinvar/clinvar_GRCh37.vcf.gz | awk '{print $5}') bytes)"

# PharmCAT (PharmGKB) for star-allele pharmacogenomics. Java 11+.
PC_VER=${PHARMCAT_VERSION:-3.4.0}
if [ ! -s tools/pharmcat/pharmcat.jar ]; then
  mkdir -p tools/pharmcat
  curl -sL -o tools/pharmcat/pharmcat.jar "https://github.com/PharmGKB/PharmCAT/releases/download/v${PC_VER}/pharmcat-${PC_VER}-all.jar"
  curl -sL -o tools/pharmcat/pharmcat_positions.vcf "https://github.com/PharmGKB/PharmCAT/releases/download/v${PC_VER}/pharmcat_positions_${PC_VER}.vcf"
fi
echo "pharmcat ok"

# Optional: ezancestry needs Python <= 3.11. If python3.11 is available, build a side venv for it.
if [ ! -x .venv-ancestry/bin/python ] && command -v python3.11 >/dev/null; then
  python3.11 -m venv .venv-ancestry && .venv-ancestry/bin/pip install -q "snps>=2.12" "git+https://github.com/arvkevi/ezancestry" "git+https://github.com/genid/Yleaf" \
    && echo "ezancestry + Yleaf ok (use --ancestry; Yleaf runs automatically for male samples)" || echo "ezancestry/Yleaf install failed (optional)"
elif [ -x .venv-ancestry/bin/python ]; then
  echo "ezancestry ok (use --ancestry)"
else
  echo "ezancestry skipped: no python3.11 on PATH (optional; brew install python@3.11)"
fi

echo
echo "Run:  .venv/bin/gnomancer /path/to/raw-file.txt [--ancestry]"
