"""The command-line tools and reference genome that sequencing analyses need, kept in $GNOMANCER_HOME.

Reading a VCF needs nothing beyond Python. Everything that touches reads (a BAM or CRAM), and Yleaf's VCF mode,
needs samtools, bcftools and friends. These are C programs distributed through bioconda, so `gnomancer setup
--wgs` builds one conda environment at ~/.gnomancer/wgs-env and gnomancer calls the binaries from there. Nothing
is installed system-wide. An existing conda, mamba or micromamba is used when there is one; otherwise the
single-file micromamba is downloaded into ~/.gnomancer/tools.

Optional callers are grouped, and each group gets its own environment (~/.gnomancer/wgs-env-<group>): their
dependencies pin different Python versions, and a caller that will not install on some platform must not take
the rest down.
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

from .paths import HOME

ENV = HOME / "wgs-env"
REF_DIR = HOME / "ref"
CHANNELS = ["-c", "conda-forge", "-c", "bioconda", "--override-channels"]

# group -> conda packages. 'core' is required for any BAM/CRAM work and for the Y call from a VCF.
GROUPS = {
    "core": ["samtools>=1.19", "bcftools>=1.19", "htslib>=1.19", "mosdepth"],
    "pgx": ["pypgx"],                               # CYP2D6 and other copy-number pharmacogenes (MIT)
    "repeats": ["expansionhunter", "stranger"],     # repeat expansions, with pathogenic thresholds
    "hla": ["t1k"],                                 # HLA typing (MIT)
    "smn": ["smncopynumbercaller", "pysam"],        # SMN1/SMN2 copy number (Illumina, non-commercial licence); the recipe omits pysam
}
# What each optional group is for, shown by `gnomancer setup --wgs`.
GROUP_LABEL = {"core": "samtools, bcftools, mosdepth", "pgx": "PyPGx (CYP2D6 and other structural pharmacogenes)",
               "repeats": "ExpansionHunter + stranger (repeat expansions)", "hla": "T1K (HLA typing)",
               "smn": "SMNCopyNumberCaller (spinal muscular atrophy carrier status; free for non-commercial use)"}

# UCSC's GRCh38 'analysis set': primary chromosomes, unplaced contigs and EBV, no alternate haplotypes. This is
# the sequence nearly every pipeline aligns to, so BAM/CRAM coordinates and contig names (chr1..chrM) match.
REFERENCE = {38: dict(url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
                      md5="6d3c82e1e12b127d526395294526b9c8", name="GRCh38.fa", size="950 MB download, 3.1 GB on disk")}


def env_dir(group: str = "core") -> Path:
    return ENV if group == "core" else HOME / f"wgs-env-{group}"


def bin_dir(group: str = "core") -> Path:
    return env_dir(group) / "bin"


def tool_bin(name: str, group: str | None = None) -> str | None:
    """Path to a sequencing tool: gnomancer's environments first (a named group, then core, then the rest), then PATH."""
    order = ([group] if group else []) + ["core"] + [g for g in GROUPS if g not in (group, "core")]
    for g in order:
        p = bin_dir(g) / name
        if p.exists():
            return str(p)
    return shutil.which(name)


def env(group: str | None = None) -> dict:
    """Environment for subprocesses: the group's tools, then the core tools, first on PATH (tools call each other by name)."""
    e = dict(os.environ)
    dirs = [str(bin_dir(g)) for g in ([group] if group and group != "core" else []) + ["core"] if bin_dir(g).exists()]
    if dirs:
        e["PATH"] = os.pathsep.join(dirs + [e.get("PATH", "")])
    return e


def have(*names: str) -> bool:
    return all(tool_bin(n) for n in names)


def conda_exe() -> str | None:
    own = HOME / "tools" / "micromamba" / "bin" / "micromamba"
    if own.exists():
        return str(own)
    for name in ("micromamba", "mamba", "conda"):
        p = shutil.which(name)
        if p:
            return p
    return None


def _fetch_micromamba() -> str:
    sysname = {"Darwin": "osx", "Linux": "linux"}.get(platform.system())
    arch = {"arm64": "arm64", "aarch64": "aarch64", "x86_64": "64", "AMD64": "64"}.get(platform.machine())
    if not sysname or not arch:
        raise RuntimeError(f"no micromamba build for {platform.system()} {platform.machine()}; install conda or mamba yourself")
    dest = HOME / "tools" / "micromamba"
    dest.mkdir(parents=True, exist_ok=True)
    url = f"https://micro.mamba.pm/api/micromamba/{sysname}-{arch}/latest"
    print(f"  micromamba: downloading {url}")
    tmp = dest / "micromamba.tar.bz2"
    urllib.request.urlretrieve(url, tmp)  # noqa: S310
    with tarfile.open(tmp, "r:bz2") as t:
        member = t.getmember("bin/micromamba")
        t.extract(member, dest, filter="data") if hasattr(tarfile, "data_filter") else t.extract(member, dest)
    tmp.unlink()
    exe = dest / "bin" / "micromamba"
    exe.chmod(0o755)
    return str(exe)


def install(groups=("core",), quiet: bool = False) -> dict:
    """Create or extend the tool environment. Returns {group: True | error text}."""
    exe = conda_exe() or _fetch_micromamba()
    out = {}
    for g in groups:
        pkgs = GROUPS[g]
        target = env_dir(g)
        verb = ["install"] if (target / "conda-meta").exists() else ["create"]
        cmd = [exe, *verb, "-y", "-p", str(target), *CHANNELS, *pkgs]
        if not quiet:
            print(f"  {g}: {GROUP_LABEL.get(g, ' '.join(pkgs))}", flush=True)
        p = subprocess.run(cmd, capture_output=True, text=True)
        out[g] = True if p.returncode == 0 else (p.stderr or p.stdout).strip()[-600:]
        if out[g] is not True and not quiet:
            print(f"  {g}: FAILED\n{out[g]}", file=sys.stderr)
    return out


def _download(url: str, dest: Path) -> None:
    """Fetch a large file. Python's bundled certificate store rejects some servers that send an incomplete
    chain (UCSC's download host is one); the system's curl validates them against the OS trust store, so it is
    the fallback. Verification is never switched off."""
    import ssl
    import urllib.error
    try:
        with urllib.request.urlopen(url) as r, open(dest, "wb") as f:  # noqa: S310
            shutil.copyfileobj(r, f, 1 << 22)
        return
    except (ssl.SSLError, urllib.error.URLError) as e:
        if not isinstance(getattr(e, "reason", e), ssl.SSLError) and not isinstance(e, ssl.SSLError):
            raise
    curl = shutil.which("curl")
    if not curl:
        raise RuntimeError(f"certificate check failed for {url} and curl is not available to retry with the system trust store")
    subprocess.run([curl, "-fL", "--retry", "3", "-sS", "-o", str(dest), url], check=True)


def reference(build: int) -> Path | None:
    """The indexed reference FASTA for a build, if it has been set up ($GNOMANCER_REF_<build> overrides)."""
    override = os.environ.get(f"GNOMANCER_REF_{build}")
    if override and Path(override).exists():
        return Path(override)
    spec = REFERENCE.get(build)
    if spec and (REF_DIR / spec["name"]).exists() and (REF_DIR / (spec["name"] + ".fai")).exists():
        return REF_DIR / spec["name"]
    return None


def fetch_reference(build: int = 38) -> Path:
    """Download, verify, unpack and index the reference. Needs samtools (for the index)."""
    import gzip
    import hashlib
    spec = REFERENCE[build]
    REF_DIR.mkdir(parents=True, exist_ok=True)
    fa = REF_DIR / spec["name"]
    if not fa.exists():
        gz = REF_DIR / (spec["name"] + ".gz.part")
        print(f"  reference GRCh{build}: downloading {spec['url']} ({spec['size']})")
        _download(spec["url"], gz)
        h = hashlib.md5()  # noqa: S324 - the publisher's checksum is MD5
        with open(gz, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != spec["md5"]:
            gz.unlink()
            raise RuntimeError(f"checksum mismatch for the GRCh{build} reference; file deleted")
        print("  reference: unpacking")
        tmp = fa.with_suffix(".fa.part")
        with gzip.open(gz, "rb") as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 22)
        tmp.rename(fa)
        gz.unlink()
    if not (REF_DIR / (spec["name"] + ".fai")).exists():
        samtools = tool_bin("samtools")
        if not samtools:
            raise RuntimeError("samtools is needed to index the reference; run `gnomancer setup --wgs` first")
        subprocess.run([samtools, "faidx", str(fa)], check=True)
    return fa


PYPGX_BUNDLE_URL = "https://github.com/sbslee/pypgx-bundle"
GENCODE_GTF = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"


def fetch_pypgx_bundle() -> Path:
    """PyPGx's reference panels and trained copy-number models (about 500 MB), at the tag matching the installed PyPGx."""
    dest = HOME / "tools" / "pypgx-bundle"
    if (dest / "cnv").exists():
        return dest
    pypgx, git = tool_bin("pypgx", "pgx"), shutil.which("git")
    if not pypgx or not git:
        raise RuntimeError("needs PyPGx and git")
    version = subprocess.run([pypgx, "-v"], capture_output=True, text=True, env=env("pgx")).stdout.split()[-1]
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  PyPGx data bundle {version}: cloning {PYPGX_BUNDLE_URL} (about 500 MB)", flush=True)
    p = subprocess.run([git, "clone", "-q", "--depth", "1", "--branch", version, PYPGX_BUNDLE_URL, str(dest)], capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip()[-300:])
    return dest


def build_t1k_index() -> Path:
    """T1K's HLA reference: allele sequences from IPD-IMGT/HLA and gene coordinates from GENCODE (GRCh38)."""
    import gzip
    dest = HOME / "tools" / "t1k"
    idx = dest / "hlaidx"
    if (idx / "hlaidx_dna_seq.fa").exists() and (idx / "hlaidx_dna_coord.fa").exists():
        return idx
    build = tool_bin("t1k-build.pl", "hla")
    if not build:
        raise RuntimeError("needs T1K")
    dest.mkdir(parents=True, exist_ok=True)
    gtf = dest / "gencode.basic.gtf"
    if not gtf.exists():
        print(f"  HLA index: downloading gene annotation {GENCODE_GTF} (30 MB)", flush=True)
        _download(GENCODE_GTF, dest / "gencode.gtf.gz")
        with gzip.open(dest / "gencode.gtf.gz", "rb") as src, open(gtf, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 22)
        (dest / "gencode.gtf.gz").unlink()
    print("  HLA index: downloading the IPD-IMGT/HLA database and building (a few minutes)", flush=True)
    for cmd in ([build, "-o", "hlaidx", "--download", "IPD-IMGT/HLA"], [build, "-o", "hlaidx", "-d", "hlaidx/hla.dat", "-g", str(gtf)]):
        p = subprocess.run(cmd, capture_output=True, text=True, cwd=str(dest), env=env("hla"))
        if p.returncode != 0:
            raise RuntimeError((p.stderr or p.stdout).strip()[-300:])
    if not (idx / "hlaidx_dna_coord.fa").exists():
        raise RuntimeError("t1k-build produced no coordinate file")
    return idx


def status() -> dict:
    return dict(env=str(ENV) if ENV.exists() else None,
                tools={n: bool(tool_bin(n)) for n in ("samtools", "bcftools", "bgzip", "mosdepth", "pypgx", "ExpansionHunter", "stranger", "run-t1k", "smn_caller.py")},
                reference={b: str(reference(b)) if reference(b) else None for b in REFERENCE})
