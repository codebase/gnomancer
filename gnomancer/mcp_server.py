"""MCP server: lets any MCP client (Claude Desktop, Claude Code, Cursor, ChatGPT…) run gnomancer and read the results.

Runs over stdio on the user's machine. Analysis happens in a subprocess so nothing the pipeline prints can touch
the protocol stream. Privacy: tools return the report's derived results and single-site lookups; no tool returns
the raw genotype file. Health and carrier tables are separate tools so an assistant only fetches them when asked.

Start:  gnomancer serve        (or: gnomancer-mcp)
Claude Code:   claude mcp add gnomancer -- gnomancer serve
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from . import __version__

try:
    from mcp.server.mcpserver import MCPServer
except ImportError as e:  # pragma: no cover
    raise SystemExit("The MCP server needs the 'mcp' package: pip install 'gnomancer[mcp]'  (or: pip install mcp)") from e

INSTRUCTIONS = """gnomancer analyzes a consumer raw DNA file locally and writes a report. Typical flow:
1. analyze(path) on the user's AncestryDNA / 23andMe / MyHeritage / FTDNA / LivingDNA / VCF file (zips are fine).
2. Read the returned at_a_glance; give the user a short recap (lineages, blood group, anything worth knowing).
3. narrative_prompt() returns a prompt; write the 250-400 word narrative it asks for and pass it to set_narrative(text)
   so it appears at the top of the user's report.
4. Answer questions with summary(), section(name), lookup(rsid). Fetch clinvar_findings() / pharmacogenomics() only
   when the user asks about health or medications; the report collapses those on purpose.
Rules: the raw file stays on the machine, but whatever these tools return is sent to your model provider, so fetch only
what the user's question needs and do not paste genotype rows anywhere. An array is probabilistic;
rare-site calls need clinical confirmation. Never diagnose. If the user has no file yet, read the
gnomancer://guide/get-your-data resource."""

server = MCPServer("gnomancer", version=__version__, instructions=INSTRUCTIONS, website_url="https://gnomancer.com")
_state: dict = {"report_dir": os.environ.get("GNOMANCER_REPORT"), "df_cache": {}}


def _report_dir(report_dir: str | None) -> Path:
    d = report_dir or _state.get("report_dir")
    if not d:
        raise ValueError("No report yet. Call analyze(path) first, or pass report_dir.")
    p = Path(d).expanduser()
    if not (p / "summary.json").exists():
        raise ValueError(f"{p} does not contain a gnomancer summary.json")
    return p


def _summary(report_dir: str | None) -> dict:
    return json.loads((_report_dir(report_dir) / "summary.json").read_text())


def _cli(args: list[str]) -> subprocess.CompletedProcess:
    exe = Path(sys.executable).parent / "gnomancer"
    cmd = [str(exe)] if exe.exists() else [sys.executable, "-m", "gnomancer"]
    return subprocess.run(cmd + args, capture_output=True, text=True)


@server.tool()
def analyze(path: str, out_dir: str | None = None, ancestry: bool = False, fast: bool = False) -> dict:
    """Run gnomancer on a raw DNA file (AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA, LivingDNA, VCF; .zip OK).
    Takes one to three minutes; the very first run also downloads ~250 MB of open-source tools and references (a few
    minutes), automatically. fast=True skips Yleaf, ClinVar and PharmCAT. Returns the output directory, report paths
    and the at-a-glance digest. gnomancer uploads nothing; the returned results go to your model provider like any tool output."""
    inp = Path(path).expanduser()
    if not inp.exists():
        raise ValueError(f"file not found: {inp}")
    from .setup_cmd import tools_present
    first_run = not any(tools_present().values())
    out = Path(out_dir).expanduser() if out_dir else inp.with_name(inp.stem.replace(".", "_") + "_analysis")
    args = [str(inp), "--out", str(out)]
    if ancestry:
        args.append("--ancestry")
    if fast:
        args += ["--no-yleaf", "--no-clinvar", "--no-pharmcat"]
    p = _cli(args)
    if p.returncode != 0:
        raise RuntimeError("gnomancer failed: " + (p.stderr or p.stdout)[-1500:])
    _state["report_dir"] = str(out)
    s = _summary(str(out))
    notes = list(s.get("notes") or [])
    if first_run:
        notes.append("First run: the open-source tools and reference files (~250 MB) were downloaded to ~/.gnomancer before analysis. Later runs are faster.")
    return dict(report_dir=str(out), report_html=str(out / "REPORT.html"), report_md=str(out / "REPORT.md"),
                source=s.get("source"), build=s.get("build"), sites=s.get("count"), sex=s.get("qc", {}).get("sex_inferred"),
                at_a_glance=s.get("glance"), notes=notes)


@server.tool()
def summary(report_dir: str | None = None) -> dict:
    """Headline results of the current (or given) report: quality, lineages, ancestry estimate, multi-site readings and
    the at-a-glance digest. No health tables; use clinvar_findings / pharmacogenomics for those."""
    s = _summary(report_dir)
    return dict(source=s.get("source"), sites=s.get("count"), qc=s.get("qc"), y=s.get("y"), mt=s.get("mt"),
                ancestry=s.get("ancestry"), readings=s.get("combos"), at_a_glance=s.get("glance"), roh_total_mb=s.get("roh_total_mb"))


@server.tool()
def section(name: str | None = None, report_dir: str | None = None) -> str:
    """Return one section of REPORT.md by (partial, case-insensitive) heading, e.g. 'Paternal', 'Marker panel',
    'PharmCAT'. With no name, lists the available section headings."""
    md = (_report_dir(report_dir) / "REPORT.md").read_text()
    parts = re.split(r"^(## .+)$", md, flags=re.M)
    heads = parts[1::2]
    if not name:
        return "\n".join(h.lstrip("# ") for h in heads)
    for h, body in zip(heads, parts[2::2]):
        if name.lower() in h.lower():
            return h + "\n" + body.strip()
    raise ValueError(f"no section matching {name!r}; headings: " + "; ".join(h.lstrip('# ') for h in heads))


@server.tool()
def lookup(rsid: str, report_dir: str | None = None) -> dict:
    """Genotype and reading for one rsID. Uses the curated panel when the site is in it; otherwise reads the single
    site from the user's raw file (nothing else is read or returned). Says so if the chip did not type the site."""
    d = _report_dir(report_dir)
    rsid = rsid.strip()
    with open(d / "marker-panel.tsv", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("rsid") == rsid:
                row = {k: r[k] for k in ("rsid", "gene", "label", "category", "genotype", "effect", "effect_count", "interpretation", "note") if k in r}
                if row.get("genotype") == "not on array":
                    row["note"] = "this vendor's chip does not type the site; no genotype exists for it"
                return row
    inp = _summary(report_dir).get("input")
    if not inp or not Path(inp).exists():
        return dict(rsid=rsid, genotype=None, note="not in the curated panel and the raw file is not available for lookup")
    df = _state["df_cache"].get(inp)
    if df is None:
        from .load import load
        df = load(inp).df
        _state["df_cache"] = {inp: df}
    if rsid in df.index:
        row = df.loc[rsid]
        g = row.genotype if isinstance(row.genotype, str) else None
        return dict(rsid=rsid, chrom=str(row.chrom), pos=int(row.pos), genotype=g or "no call",
                    note="raw site, not in the curated panel: no interpretation attached; alleles are GRCh37 forward strand")
    return dict(rsid=rsid, genotype=None, note="not typed on this chip")


@server.tool()
def clinvar_findings(report_dir: str | None = None) -> dict:
    """ClinVar scan results (pathogenic / likely-pathogenic alleles carried, with plausibility confidence and flags).
    Health-sensitive: call only when the user asks about disease risk or carrier status."""
    s = _summary(report_dir)
    return s.get("clinvar") or {}


@server.tool()
def pharmacogenomics(report_dir: str | None = None) -> dict:
    """PharmCAT star-allele diplotypes and phenotypes per gene, with coverage. Provisional: untyped positions are
    assumed reference. Call when the user asks about medications or drug metabolism."""
    s = _summary(report_dir)
    return s.get("pharmcat") or {}


@server.tool()
def clinvar_watch(report_dir: str | None = None, force: bool = False) -> str:
    """Re-check ClinVar for this report: downloads a new ClinVar release if NCBI published one, rescans the user's
    sites and reports only classification changes (newly pathogenic, downgraded, review-status changes). The first
    run writes a baseline. Health-sensitive: run when the user asks."""
    d = _report_dir(report_dir)
    args = ["watch", str(d)] + (["--force"] if force else [])
    p = _cli(args)
    if p.returncode != 0:
        raise RuntimeError((p.stderr or p.stdout)[-1200:])
    return (p.stdout or "").strip() or (p.stderr or "").strip()


@server.tool()
def family(paths: list[str], names: list[str] | None = None, out_dir: str | None = None, haplogroups: bool = True) -> dict:
    """Compare two or more relatives' raw DNA files: shared DNA per pair with a relationship estimate, shared maternal/
    paternal lines, side-by-side panel genotypes and, for parent-parent-child sets, which parent each allele came from
    plus a Mendelian-consistency check. Runs locally. Returns the pair table and report paths."""
    ps = [Path(p).expanduser() for p in paths]
    for p in ps:
        if not p.exists():
            raise ValueError(f"file not found: {p}")
    out = Path(out_dir).expanduser() if out_dir else ps[0].with_name("family_analysis")
    args = ["family", *map(str, ps), "--out", str(out)] + (["--names", ",".join(names)] if names else []) + ([] if haplogroups else ["--no-haplogroups"])
    r = _cli(args)
    if r.returncode != 0:
        raise RuntimeError("gnomancer family failed: " + (r.stderr or r.stdout)[-1500:])
    md = (out / "family-report.md").read_text()
    table = md.split("## Who shares what", 1)[1].split("\n\n", 2)[1] if "## Who shares what" in md else ""
    return dict(report_html=str(out / "family-report.html"), report_md=str(out / "family-report.md"), pairs=table)


@server.tool()
def narrative_prompt(report_dir: str | None = None) -> str:
    """The prompt for writing the plain-language narrative that goes at the top of the report. Follow it exactly,
    then pass your text to set_narrative."""
    return (_report_dir(report_dir) / "summary-prompt.md").read_text()


@server.tool()
def set_narrative(text: str, report_dir: str | None = None) -> dict:
    """Embed a narrative summary (yours, written from narrative_prompt) at the top of the report. Saves summary.md
    and regenerates REPORT.html / REPORT.md with the same options as the original run."""
    d = _report_dir(report_dir)
    (d / "summary.md").write_text(text.strip() + "\n")
    s = _summary(str(d))
    args = list(s.get("args") or [])
    if not args:
        raise RuntimeError("this report predates argument recording; re-run gnomancer with --summary summary.md")
    kept, skip = [], False
    for a_ in args:  # drop any earlier --summary FILE and --ai so the rerun only embeds this text
        if skip:
            skip = False
            continue
        if a_ == "--summary":
            skip = True
            continue
        if a_.startswith("--summary=") or a_ == "--ai":
            continue
        kept.append(a_)
    p = _cli(kept + ["--summary", str(d / "summary.md")])
    if p.returncode != 0:
        raise RuntimeError("regeneration failed: " + (p.stderr or p.stdout)[-1500:])
    return dict(report_html=str(d / "REPORT.html"), embedded=True)


@server.resource("gnomancer://guide/get-your-data", name="How to download raw DNA data", mime_type="text/markdown")
def guide_get_your_data() -> str:
    """Per-vendor steps for downloading the raw genotype file."""
    return """# Getting your raw DNA data (summary; full steps with official links at https://gnomancer.com/get-your-data/)
- AncestryDNA: DNA results > Settings (gear) > Download DNA data > confirm password > click the emailed link (expires in 7 days). Beware the adjacent Delete option.
- 23andMe: name menu > Settings > 23andMe Data > View > Download Raw Data (plain raw genotype). Emailed when ready; can take days to weeks.
- MyHeritage: DNA > Manage DNA kits > three dots > Download kit. Link lasts 24 h; computer or Android only; kit manager only; uploaded kits cannot be re-downloaded.
- FamilyTreeDNA: Results & Tools > Autosomal DNA > Download Raw Data > Build 37 Concatenated (.csv.gz).
- LivingDNA: left sidebar > Download Raw Data > consent box > autosomal (family ancestry) file.
Do it on a computer. Keep the file private. gnomancer reads the zip directly and detects the vendor."""


@server.resource("gnomancer://report/markdown", name="Current report (Markdown)", mime_type="text/markdown")
def report_markdown() -> str:
    """The full REPORT.md of the current report."""
    return (_report_dir(None) / "REPORT.md").read_text()


def main() -> None:
    server.run("stdio")


if __name__ == "__main__":
    main()
