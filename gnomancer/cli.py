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
from .report_html import write_html

REPO = Path(__file__).resolve().parent.parent


def run_ezancestry(inp: Path, loaded) -> dict | None:
    """Run ezancestry in a separate interpreter (it needs Python <=3.11) if .venv-ancestry exists, else in-process."""
    import json as _json
    import subprocess
    code = ("import json,sys\nfrom snps import SNPs\nr=SNPs(sys.argv[1]).predict_ancestry()\n"
            "print(json.dumps({k:v for k,v in r.items() if k!='ezancestry_df'}, default=str))")
    venv_py = REPO / ".venv-ancestry" / "bin" / "python"
    py = str(venv_py) if venv_py.exists() else sys.executable
    try:
        p = subprocess.run([py, "-c", code, str(inp)], capture_output=True, text=True, timeout=1800)
        if p.returncode != 0:
            raise RuntimeError(p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "unknown error")
        return _json.loads(p.stdout.strip().splitlines()[-1])
    except Exception as e:  # noqa: BLE001
        loaded.notes.append(f"ezancestry failed: {e}")
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gnomancer", description="Analyze a consumer raw DNA file locally.")
    ap.add_argument("input", help="raw data file (.txt/.csv/.tsv/.vcf, optionally .zip/.gz)")
    ap.add_argument("--out", default=None, help="output directory (default: <input>_analysis)")
    ap.add_argument("--markers", default=str(REPO / "data" / "markers.tsv"), help="marker panel TSV")
    ap.add_argument("--no-y", action="store_true", help="skip Y haplogroup")
    ap.add_argument("--no-mt", action="store_true", help="skip mtDNA haplogroup")
    ap.add_argument("--ancestry", action="store_true", help="run ezancestry population prediction (uses .venv-ancestry if present)")
    ap.add_argument("--no-clinvar", action="store_true", help="skip the ClinVar scan")
    ap.add_argument("--no-yleaf", action="store_true", help="skip the Yleaf second-opinion Y call")
    ap.add_argument("--no-pharmcat", action="store_true", help="skip PharmCAT")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args(argv)

    inp = Path(a.input)
    out = Path(a.out) if a.out else inp.with_name(inp.stem.replace(".", "_") + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

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
    cv_text = clinvar_mod.interpret(cv) if not a.no_clinvar else cv["note"]
    pc = dict(available=False, calls=None, error="skipped (--no-pharmcat)", stats={})
    if not a.no_pharmcat:
        print("[7/8] PharmCAT pharmacogenomics ...", file=sys.stderr)
        pc = pharmcat_mod.run(df, out)
        if pc.get("calls") is not None and not pc["calls"].empty:
            pc["calls"].to_csv(out / "pharmcat-calls.tsv", sep="\t", index=False)
    pc_text = pharmcat_mod.interpret(pc) if not a.no_pharmcat else pc["error"]
    ancestry = None
    if a.ancestry:
        print("[8/8] ancestry prediction (ezancestry) ...", file=sys.stderr)
        ancestry = run_ezancestry(inp, loaded)

    extra = dict(clinvar=cv, clinvar_text=cv_text, pharmcat=pc, pharmcat_text=pc_text)
    rep = write(out / "REPORT.md", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, extra)
    html = write_html(out / "REPORT.html", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, input_name=inp.name, extra=extra)
    json.dump(dict(input=str(inp), source=loaded.source, build=loaded.build, count=loaded.count, qc=qc,
                   roh_total_mb=float(roh.length_mb.sum()) if not roh.empty else 0.0, y=y, mt=mt,
                   combos=combo, ancestry=ancestry, notes=loaded.notes,
                   clinvar=dict(sites_in_clinvar=cv.get("sites_in_clinvar"), pathogenic_sites_typed=cv.get("pathogenic_sites_typed"),
                                hits=(cv["hits"].to_dict("records") if cv.get("hits") is not None and not cv["hits"].empty else [])),
                   pharmcat=dict(stats={k: v for k, v in pc.get("stats", {}).items() if k != "genes"}, error=pc.get("error"),
                                 calls=(pc["calls"].to_dict("records") if pc.get("calls") is not None and not pc["calls"].empty else []))),
              open(out / "summary.json", "w"), indent=1, default=str)
    print(f"\nReport: {html}\n        {rep}", file=sys.stderr)
    if y.get("haplogroup"):
        print(f"Y: {y['haplogroup']}   ", end="", file=sys.stderr)
    if mt.get("haplogroup"):
        print(f"MT: {mt['haplogroup']}", end="", file=sys.stderr)
    print(file=sys.stderr)
