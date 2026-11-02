"""Command-line entry point: python -m gnomancer <raw file> [--out DIR]"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .load import load
from .qc import chromosome_table, summary
from .roh import find_roh, interpret
from .panel import load_markers, evaluate, combos
from .haplogroups import call_y, call_mt
from .report import write
from . import clinvar as clinvar_mod
from . import pharmcat as pharmcat_mod
from . import yleaf as yleaf_mod
from . import summary as summary_mod
from . import pgxcard
from .report_html import write_html

from .paths import DATA, ancestry_python


def run_ezancestry(inp: Path, loaded) -> dict | None:
    """Run ezancestry in a separate interpreter (it needs Python <=3.11) if .venv-ancestry exists, else in-process."""
    import json as _json
    import subprocess
    code = ("import json,sys\nfrom snps import SNPs\nr=SNPs(sys.argv[1]).predict_ancestry()\n"
            "print(json.dumps({k:v for k,v in r.items() if k!='ezancestry_df'}, default=str))")
    venv_py = ancestry_python()
    py = str(venv_py) if venv_py else sys.executable
    try:
        p = subprocess.run([py, "-c", code, str(inp)], capture_output=True, text=True, timeout=1800)
        if p.returncode != 0:
            raise RuntimeError(p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "unknown error")
        return _json.loads(p.stdout.strip().splitlines()[-1])
    except Exception as e:  # noqa: BLE001
        loaded.notes.append(f"ezancestry failed: {e}")
        return None


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["serve"]:
        from .mcp_server import main as serve
        return serve()
    if argv[:1] == ["setup"]:
        from .setup_cmd import run as setup
        return setup(argv[1:])
    if argv[:1] == ["watch"]:
        from .watch import run as watch
        return watch(argv[1:])
    ap = argparse.ArgumentParser(prog="gnomancer", description="Analyze a consumer raw DNA file locally. Subcommands: `gnomancer setup` (download tools), `gnomancer serve` (MCP server), `gnomancer watch REPORT_DIR` (ClinVar reclassification check).")
    ap.add_argument("input", help="raw data file (.txt/.csv/.tsv/.vcf, optionally .zip/.gz)")
    ap.add_argument("--out", default=None, help="output directory (default: <input>_analysis)")
    ap.add_argument("--markers", default=str(DATA / "markers.tsv"), help="marker panel TSV")
    ap.add_argument("--no-y", action="store_true", help="skip Y haplogroup")
    ap.add_argument("--no-mt", action="store_true", help="skip mtDNA haplogroup")
    ap.add_argument("--ancestry", action="store_true", help="run ezancestry population prediction (uses .venv-ancestry if present)")
    ap.add_argument("--no-clinvar", action="store_true", help="skip the ClinVar scan")
    ap.add_argument("--no-setup", action="store_true", help="do not auto-download tools on first run")
    ap.add_argument("--no-yleaf", action="store_true", help="skip the Yleaf second-opinion Y call")
    ap.add_argument("--no-pharmcat", action="store_true", help="skip PharmCAT")
    ap.add_argument("--summary", default=None, help="Markdown/text file with a narrative summary (e.g. written by your AI agent) to show at the top of the report")
    ap.add_argument("--ai", action="store_true", help="have Claude write the narrative summary (sends the structured results, not the raw file, to the Anthropic API; needs credentials)")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args(argv)
    run_args = argv

    inp = Path(a.input)
    out = Path(a.out) if a.out else inp.with_name(inp.stem.replace(".", "_") + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

    if not a.no_setup:
        from .setup_cmd import ensure_setup
        ensure_setup()
    print(f"[1/6] loading {inp} ...", file=sys.stderr)
    loaded = load(str(inp))
    df = loaded.df
    print(f"      {loaded.source}, build {loaded.build}, {loaded.count:,} SNPs, sex {loaded.sex or '?'}", file=sys.stderr)

    print("[2/6] quality metrics ...", file=sys.stderr)
    ct = chromosome_table(df)
    qc = summary(df, ct, loaded.sex)
    ct.to_csv(out / "chromosome-qc.tsv", sep="\t", index=False)

    print("[3/6] runs of homozygosity ...", file=sys.stderr)
    roh = find_roh(df)
    roh.to_csv(out / "roh.tsv", sep="\t", index=False)
    roh_text = interpret(roh)

    print("[4/6] marker panel ...", file=sys.stderr)
    markers = load_markers(Path(a.markers))
    panel = evaluate(df, markers)
    panel.to_csv(out / "marker-panel.tsv", sep="\t", index=False)
    combo = combos(df)

    is_male = qc["sex_inferred"] == "Male"
    y = dict(error="skipped (sample is not male)" if not is_male else "skipped (--no-y)")
    if is_male and not a.no_y:
        print("[5/6] Y haplogroup (yhaplo) ...", file=sys.stderr)
        y = call_y(df, out)
        if not a.no_yleaf:
            print("      second opinion on current trees (Yleaf) ...", file=sys.stderr)
            y["yleaf"] = yleaf_mod.run(df, out)
    mt = dict(error="skipped (--no-mt)")
    if not a.no_mt:
        print("[5/6] mtDNA haplogroup (haplogrep3) ...", file=sys.stderr)
        mt = call_mt(df, out, loaded.build_original)

    cv = dict(available=False, hits=None, note="skipped (--no-clinvar)")
    if not a.no_clinvar:
        print("[6/8] ClinVar scan of every typed site ...", file=sys.stderr)
        cv = clinvar_mod.scan(df)
        if cv["available"] and cv["hits"] is not None and not cv["hits"].empty:
            cv["hits"].to_csv(out / "clinvar-hits.tsv", sep="\t", index=False)
        if cv["available"]:
            from .watch import write_snapshot
            write_snapshot(out, cv, str(inp.resolve()))
    cv_text = clinvar_mod.interpret(cv) if not a.no_clinvar else cv["note"]
    pc = dict(available=False, calls=None, error="skipped (--no-pharmcat)", stats={})
    if not a.no_pharmcat:
        print("[7/8] PharmCAT pharmacogenomics ...", file=sys.stderr)
        pc = pharmcat_mod.run(df, out)
        if pc.get("calls") is not None and not pc["calls"].empty:
            pc["calls"].to_csv(out / "pharmcat-calls.tsv", sep="\t", index=False)
    pc_text = pharmcat_mod.interpret(pc) if not a.no_pharmcat else pc["error"]
    card = pgxcard.write(out, pc, panel, loaded.source, inp.name)
    ancestry = None
    if a.ancestry:
        print("[8/8] ancestry prediction (ezancestry) ...", file=sys.stderr)
        ancestry = run_ezancestry(inp, loaded)

    glance = summary_mod.at_a_glance(loaded, qc, roh, roh_text, combo, y, mt, ancestry, panel, cv, pc)
    prompt = summary_mod.build_prompt(glance, combo, y, mt, ancestry, panel, cv, pc, qc)
    summary_mod.write_prompt_file(out, prompt)
    narrative = None
    if a.summary:
        narrative = Path(a.summary).read_text().strip()
    elif a.ai:
        print("      narrative summary (Claude) ...", file=sys.stderr)
        try:
            narrative = summary_mod.ai_summary(prompt)
            (out / "summary.md").write_text(narrative + "\n")
        except Exception as e:  # noqa: BLE001
            loaded.notes.append(f"AI summary failed: {e}")
    pop = (ancestry or {}).get("superpopulation_code") or "EUR"
    extra = dict(clinvar=cv, clinvar_text=cv_text, pharmcat=pc, pharmcat_text=pc_text, glance=glance, narrative=narrative, pop=pop, card=card)
    rep = write(out / "REPORT.md", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, extra)
    html = write_html(out / "REPORT.html", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, input_name=inp.name, extra=extra)
    json.dump(dict(input=str(inp.resolve()), args=[str(Path(x).resolve()) if Path(x).exists() else x for x in run_args],
                   report_html=str((out / "REPORT.html").resolve()), report_md=str((out / "REPORT.md").resolve()), glance=glance,
                   source=loaded.source, build=loaded.build, count=loaded.count, qc=qc,
                   roh_total_mb=float(roh.length_mb.sum()) if not roh.empty else 0.0, y=y, mt=mt,
                   combos=combo, ancestry=ancestry, notes=loaded.notes,
                   clinvar=dict(sites_in_clinvar=cv.get("sites_in_clinvar"), pathogenic_sites_typed=cv.get("pathogenic_sites_typed"),
                                hits=(cv["hits"].to_dict("records") if cv.get("hits") is not None and not cv["hits"].empty else [])),
                   pharmcat=dict(stats={k: v for k, v in pc.get("stats", {}).items() if k != "genes"}, error=pc.get("error"),
                                 calls=(pc["calls"].to_dict("records") if pc.get("calls") is not None and not pc["calls"].empty else []))),
              open(out / "summary.json", "w"), indent=1, default=str)
    print(f"\nReport: {html}\n        {rep}" + (f"\nPGx card: {card['pdf'] or card['html']}" if card.get("rows") else ""), file=sys.stderr)
    if y.get("haplogroup"):
        print(f"Y: {y['haplogroup']}   ", end="", file=sys.stderr)
    if mt.get("haplogroup"):
        print(f"MT: {mt['haplogroup']}", end="", file=sys.stderr)
    print(file=sys.stderr)
