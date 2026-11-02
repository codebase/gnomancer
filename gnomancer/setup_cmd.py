"""`gnomancer setup`: download the open-source tools and reference files into $GNOMANCER_HOME (default ~/.gnomancer)."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

from .paths import HOME

HAPLOGREP_VERSION = "3.3.2"
PHARMCAT_VERSION = "3.4.0"
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
YHAPLO_SRC = "git+https://github.com/23andMe/yhaplo"
EZANCESTRY_SRC = "git+https://github.com/arvkevi/ezancestry"
YLEAF_SRC = "git+https://github.com/genid/Yleaf"


def _dl(url: str, dest: Path, label: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  {label}: already present")
        return
    print(f"  {label}: downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:  # noqa: S310
        shutil.copyfileobj(r, f)
    tmp.rename(dest)


def _java_ok() -> str | None:
    for jh in (os.environ.get("JAVA_HOME"), "/opt/homebrew/opt/openjdk", "/usr/local/opt/openjdk", "/usr/lib/jvm/default-java"):
        if jh and Path(jh, "bin", "java").exists():
            return str(Path(jh, "bin", "java"))
    return shutil.which("java")


def run(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gnomancer setup", description="Download tools and reference files (once).")
    ap.add_argument("--home", default=str(HOME), help="where to put them (default: $GNOMANCER_HOME or ~/.gnomancer)")
    ap.add_argument("--no-clinvar", action="store_true", help="skip the ~200 MB ClinVar download")
    ap.add_argument("--no-ancestry", action="store_true", help="skip the optional Python 3.11 venv for ezancestry and Yleaf")
    a = ap.parse_args(argv)
    home = Path(a.home).expanduser()
    tools = home / "tools"
    tools.mkdir(parents=True, exist_ok=True)
    print(f"gnomancer setup -> {home}")

    print("yhaplo (Y haplogroups)")
    try:
        import yhaplo  # noqa: F401
        print("  already installed")
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", YHAPLO_SRC], check=False)

    print("Haplogrep3 (mtDNA haplogroups)")
    hg = tools / "haplogrep3"
    if not (hg / "haplogrep3").exists():
        tgz = tools / "haplogrep3.tar.gz"
        _dl(f"https://github.com/genepi/haplogrep3/releases/download/v{HAPLOGREP_VERSION}/haplogrep3-{HAPLOGREP_VERSION}-linux.tar.gz", tgz, "haplogrep3")
        hg.mkdir(exist_ok=True)
        with tarfile.open(tgz) as t:
            t.extractall(hg, filter="data")
        tgz.unlink()
        (hg / "haplogrep3").chmod(0o755)
    else:
        print("  already present")

    print("PharmCAT (pharmacogenomics)")
    _dl(f"https://github.com/PharmGKB/PharmCAT/releases/download/v{PHARMCAT_VERSION}/pharmcat-{PHARMCAT_VERSION}-all.jar", tools / "pharmcat" / "pharmcat.jar", "pharmcat.jar")
    _dl(f"https://github.com/PharmGKB/PharmCAT/releases/download/v{PHARMCAT_VERSION}/pharmcat_positions_{PHARMCAT_VERSION}.vcf", tools / "pharmcat" / "pharmcat_positions.vcf", "positions")

    if not a.no_clinvar:
        print("ClinVar (GRCh37, ~200 MB)")
        _dl(CLINVAR_URL, tools / "clinvar" / "clinvar_GRCh37.vcf.gz", "clinvar")

    java = _java_ok()
    print(f"Java: {java or 'NOT FOUND (Haplogrep3 and PharmCAT will be skipped; install Java 11+, e.g. brew install openjdk)'}")

    if not a.no_ancestry:
        print("Optional: ezancestry + Yleaf (need Python 3.11)")
        venv = home / "venv-ancestry"
        py311 = shutil.which("python3.11")
        if (venv / "bin" / "python").exists():
            print("  already present")
        elif py311:
            subprocess.run([py311, "-m", "venv", str(venv)], check=True)
            r = subprocess.run([str(venv / "bin" / "pip"), "install", "-q", "snps>=2.12", EZANCESTRY_SRC, YLEAF_SRC], capture_output=True, text=True)
            print("  ok (use --ancestry; Yleaf runs for male samples)" if r.returncode == 0 else "  failed (optional): " + r.stderr[-300:])
        else:
            print("  skipped: python3.11 not on PATH (brew install python@3.11 to enable)")
    print("done. Run:  gnomancer /path/to/raw-file.zip")
    return 0
