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
from .haplogroups import call_y, call_mt, call_mt_seq
from .report import write
from . import clinvar as clinvar_mod
from . import pharmcat as pharmcat_mod
from . import yleaf as yleaf_mod
from . import summary as summary_mod
from . import pgxcard
from .report_html import write_html

from .paths import DATA, ancestry_python


def _stepper():
    """Numbered progress labels; how many steps run depends on sex and flags, so no denominator is promised."""
    n = iter(range(1, 20))
    return lambda: str(next(n))


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
    if argv[:1] == ["family"]:
        from .family import run as family
        return family(argv[1:])
    ap = argparse.ArgumentParser(prog="gnomancer", description="Analyze a consumer raw DNA file locally. Subcommands: `gnomancer setup` (download tools), `gnomancer serve` (MCP server), `gnomancer watch REPORT_DIR` (ClinVar reclassification check), `gnomancer family KIT1 KIT2 ...` (compare relatives).")
    ap.add_argument("input", help="raw data file (.txt/.csv/.tsv/.vcf, optionally .zip/.gz)")
    ap.add_argument("--out", default=None, help="output directory (default: <input>_analysis)")
    ap.add_argument("--markers", default=str(DATA / "markers.tsv"), help="marker panel TSV")
    ap.add_argument("--no-y", action="store_true", help="skip Y haplogroup")
    ap.add_argument("--no-mt", action="store_true", help="skip mtDNA haplogroup")
    ap.add_argument("--ancestry", action="store_true", help="run ezancestry population prediction (uses .venv-ancestry if present)")
    ap.add_argument("--no-clinvar", action="store_true", help="skip the ClinVar scan")
    ap.add_argument("--no-ancestry", action="store_true", help="sequencing input: skip the population estimate that runs automatically when ezancestry is installed")
    ap.add_argument("--no-pgs", action="store_true", help="sequencing input: skip polygenic scores")
    ap.add_argument("--no-setup", action="store_true", help="do not auto-download tools on first run")
    ap.add_argument("--no-yleaf", action="store_true", help="skip the Yleaf second-opinion Y call")
    ap.add_argument("--no-pharmcat", action="store_true", help="skip PharmCAT")
    ap.add_argument("--summary", default=None, help="Markdown/text file with a narrative summary (e.g. written by your AI agent) to show at the top of the report")
    ap.add_argument("--ai", action="store_true", help="have Claude write the narrative summary (sends the structured results, not the raw file, to the Anthropic API; needs credentials)")
    ap.add_argument("--sv", action="append", default=[], metavar="VCF", help="with a sequencing VCF: the vendor's structural-variant VCF (repeatable)")
    ap.add_argument("--cnv", action="append", default=[], metavar="VCF", help="with a sequencing VCF: the vendor's copy-number VCF (repeatable)")
    ap.add_argument("--bam", default=None, metavar="BAM_OR_CRAM", help="with a sequencing VCF: the aligned reads, for coverage, confirming reference "
                    "positions, the Y lineage and mitochondrial heteroplasmy (needs `gnomancer setup --wgs`)")
    ap.add_argument("--mito", default=None, metavar="VCF", help="with a sequencing VCF: the vendor's mitochondrial heteroplasmy VCF, if you have it and no BAM")
    ap.add_argument("--callers", default="all", metavar="LIST", help="with --bam: which read-level callers to run, comma-separated from "
                    "repeats,cyp2d6,hla,smn (default all that are installed; 'none' to skip). Together they add about five minutes on a 30x genome")
    ap.add_argument("--build", type=int, choices=(37, 38), default=None, help="genome build of a sequencing VCF, if its header does not say")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args(argv)
    run_args = argv

    inp = Path(a.input)
    out = Path(a.out) if a.out else inp.with_name(inp.stem.replace(".", "_") + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

    if not a.no_setup:
        from .setup_cmd import ensure_setup
        ensure_setup()
    _step = _stepper()
    print(f"[{_step()}] loading {inp} ...", file=sys.stderr)
    from . import seq as seq_mod
    sn = seq_mod.detect(str(inp))
    if sn and not a.no_clinvar and not a.no_setup and (a.build or sn.build) in (37, 38):
        from .setup_cmd import ensure_clinvar
        try:
            ensure_clinvar(a.build or sn.build)
        except Exception as e:  # noqa: BLE001
            print(f"      could not fetch ClinVar for this build ({e}); the ClinVar scan will be skipped", file=sys.stderr)
    if sn:
        print(f"      sequencing VCF ({sn.why}); reading every record, this takes a minute or two ...", file=sys.stderr)
        loaded = seq_mod.load_seq(str(inp), sn, build=a.build, use_clinvar=not a.no_clinvar, use_pharmcat=not a.no_pharmcat)
    else:
        loaded = load(str(inp))
    df, seq = loaded.df, loaded.seq
    print(f"      {loaded.source}, build {loaded.build}, {loaded.count:,} {'variants' if seq else 'SNPs'}, sex {loaded.sex or '?'}", file=sys.stderr)

    reads = None
    if a.bam and not seq:
        loaded.notes.append("--bam is read only alongside a sequencing VCF; it was ignored.")
    elif a.bam:
        from . import bam as bam_mod
        print(f"[{_step()}] reads: coverage, and every named position the VCF is silent about ...", file=sys.stderr)
        try:
            reads = bam_mod.open_bam(a.bam)
            if reads.build != loaded.build:
                raise bam_mod.BamError(f"{reads.path.name} is aligned to GRCh{reads.build} but the VCF is on GRCh{loaded.build}")
            reads.cov = bam_mod.coverage(reads)
            reads.confirm = bam_mod.confirm_reference(reads, seq, seq.named)
            seq_mod.refresh(loaded)
            df = loaded.df
            loaded.notes += reads.notes
        except bam_mod.BamError as e:
            loaded.notes.append(f"The BAM/CRAM could not be used: {e}")
            reads = None

    print(f"[{_step()}] quality metrics ...", file=sys.stderr)
    if seq:
        from . import seqqc
        ct, qc = seqqc.chromosome_table(seq), seqqc.summary(seq)
        if reads:
            qc["reads"] = dict(file=reads.path.name, **{k: v for k, v in reads.cov.items() if k != "per_chrom"}, **{f"sites_{k}": v for k, v in reads.confirm.items()})
            qc["flags"] += bam_mod.flags(reads, seq)
    else:
        ct = chromosome_table(df)
        qc = summary(df, ct, loaded.sex)
    ct.to_csv(out / "chromosome-qc.tsv", sep="\t", index=False)

    print(f"[{_step()}] runs of homozygosity ...", file=sys.stderr)
    if seq and not qc["whole_genome"]:
        roh, roh_text = find_roh(df.iloc[0:0]), "Not assessed: runs of homozygosity need a whole genome, and this file covers only part of one."
    else:
        roh = seqqc.find_roh(seq) if seq else find_roh(df)
        roh_text = interpret(roh)
    roh.to_csv(out / "roh.tsv", sep="\t", index=False)

    print(f"[{_step()}] marker panel ...", file=sys.stderr)
    markers = load_markers(Path(a.markers))
    panel = evaluate(df, markers)
    panel.to_csv(out / "marker-panel.tsv", sep="\t", index=False)
    combo = combos(df)

    is_male = qc["sex_inferred"] == "Male"
    y = dict(error="skipped (sample is not male)" if not is_male else "skipped (--no-y)")
    if is_male and not a.no_y and seq:
        print(f"[{_step()}] Y haplogroup (Yleaf{', from the reads' if reads else ''}) ...", file=sys.stderr)
        y = bam_mod.yleaf(reads, out) if reads else yleaf_mod.run_seq(seq, out, whole_genome=qc["whole_genome"])
        if reads and y.get("error"):          # fall back to the VCF if the read-based call failed
            first = y["error"]
            y = yleaf_mod.run_seq(seq, out, whole_genome=qc["whole_genome"])
            loaded.notes.append(f"Yleaf could not use the BAM/CRAM ({first.strip().splitlines()[-1][:200]}); the Y call is from the VCF.")
    elif is_male and not a.no_y:
        print(f"[{_step()}] Y haplogroup (yhaplo) ...", file=sys.stderr)
        y = call_y(df, out)
        if not a.no_yleaf:
            print("      second opinion on current trees (Yleaf) ...", file=sys.stderr)
            y["yleaf"] = yleaf_mod.run(df, out)
    mt = dict(error="skipped (--no-mt)")
    if not a.no_mt:
        print(f"[{_step()}] mtDNA haplogroup (haplogrep3) ...", file=sys.stderr)
        if seq and a.mito:
            try:
                added = seq_mod.merge_mito_vcf(seq, a.mito)
                loaded.notes.append(f"{added} mitochondrial variant(s) from {Path(a.mito).name} were added to the calls in the main VCF.")
            except Exception as e:  # noqa: BLE001
                loaded.notes.append(f"The mitochondrial VCF could not be used: {e}")
        if reads:
            het = bam_mod.mito_heteroplasmy(reads)
            bam_mod.merge_mito(seq, het)
        mt = call_mt_seq(seq, out) if seq else call_mt(df, out, loaded.build_original)
        if reads and seq:
            mt["read_level"] = dict(depth=het.get("depth"), error=het.get("error"), floor=bam_mod.MT_MIN_SHARE)

    cv = dict(available=False, hits=None, note="skipped (--no-clinvar)")
    if not a.no_clinvar:
        print(f"[{_step()}] ClinVar scan of every {'variant' if seq else 'typed site'} ...", file=sys.stderr)
        cv = clinvar_mod.scan_seq(seq) if seq else clinvar_mod.scan(df)
        if cv["available"] and cv["hits"] is not None and not cv["hits"].empty:
            cv["hits"].to_csv(out / "clinvar-hits.tsv", sep="\t", index=False)
        if cv["available"]:
            from .watch import write_snapshot
            write_snapshot(out, cv, str(inp.resolve()))
    cv_text = cv["note"] if a.no_clinvar else clinvar_mod.interpret_seq(cv, qc["whole_genome"]) if seq else clinvar_mod.interpret(cv)
    called = {}
    if reads and a.callers != "none":
        import time
        from . import callers as callers_mod
        wanted = None if a.callers == "all" else {x.strip() for x in a.callers.split(",")}
        for name, fn, _ in callers_mod.ALL:
            if wanted is not None and name not in wanted:
                continue
            print(f"[{_step()}] {dict(repeats='repeat expansions (ExpansionHunter)', cyp2d6='CYP2D6 from read depth (PyPGx)', hla='HLA typing (T1K)', smn='SMN1/SMN2 copy number')[name]} ...", file=sys.stderr)
            t0 = time.time()
            try:
                called[name] = fn(reads, out, qc["sex_inferred"]) if name == "repeats" else fn(reads, out)
            except Exception as e:  # noqa: BLE001 - one failing caller must not lose the rest of the report
                called[name] = dict(available=False, error=f"{type(e).__name__}: {e}", section=None, what=name)
            called[name]["minutes"] = round((time.time() - t0) / 60, 1)
            if not called[name]["available"]:
                print(f"      not run: {called[name]['error']}", file=sys.stderr)
    outside = {c["gene"]: c["genotype"] for c in (called.get("cyp2d6") or {}).get("calls", []) if c.get("genotype") and "Indeterminate" not in c["genotype"]}
    pc = dict(available=False, calls=None, error="skipped (--no-pharmcat)", stats={})
    if not a.no_pharmcat:
        print(f"[{_step()}] PharmCAT pharmacogenomics ...", file=sys.stderr)
        pc = pharmcat_mod.run(df, out, seq=seq, whole_genome=bool(seq and qc["whole_genome"]), outside=outside or None)
        if pc.get("calls") is not None and not pc["calls"].empty:
            pc["calls"].to_csv(out / "pharmcat-calls.tsv", sep="\t", index=False)
    pc_text = pc["error"] if a.no_pharmcat else pharmcat_mod.interpret_seq(pc, seq) if seq else pharmcat_mod.interpret(pc)
    card = pgxcard.write(out, pc, panel, loaded.source, inp.name)
    ancestry = None
    if a.ancestry or (seq and not a.no_ancestry and ancestry_python()):
        print(f"[{_step()}] ancestry prediction (ezancestry) ...", file=sys.stderr)
        ancestry = run_ezancestry(seq_mod.ancestry_file(loaded, out) if seq else inp, loaded)

    sections, more = [], {}
    if a.sv or a.cnv:
        if not seq:
            loaded.notes.append("--sv and --cnv are read only alongside a sequencing VCF; they were ignored.")
        else:
            print(f"[{_step()}] structural and copy-number variants ...", file=sys.stderr)
            from . import sv as sv_mod
            more["sv"] = sv_mod.analyze(a.sv + a.cnv, loaded.build, qc["sex_inferred"], out)
            sections.append(more["sv"]["section"])

    for name, _, _ in (callers_mod.ALL if called else ()):
        r = called.get(name)
        if r and r.get("section"):
            sections.append(r["section"])
        elif r:
            loaded.notes.append(f"{name} was not run: {r['error']}")
    if seq and not a.no_pgs:
        print(f"[{_step()}] polygenic scores ...", file=sys.stderr)
        from . import pgs as pgs_mod
        more["pgs"] = pgs_mod.score(seq, qc["whole_genome"], qc["sex_inferred"], (ancestry or {}).get("superpopulation_code"))
        if more["pgs"]["section"]:
            sections.append(more["pgs"]["section"])

    glance = summary_mod.at_a_glance(loaded, qc, roh, roh_text, combo, y, mt, ancestry, panel, cv, pc)
    if "sv" in more:
        notable, reassuring = sv_mod.glance_bullets(more["sv"])
        glance["notable"] += notable
        glance["reassuring"] += reassuring
    for name, _, bullets in (callers_mod.ALL if called else ()):
        if name == "cyp2d6" and "CYP2D6" in (pc.get("outside") or {}):
            continue                      # PharmCAT took the diplotype and already reports it with its phenotype
        if called.get(name):
            notable, reassuring = bullets(called[name])
            glance["notable"] += notable
            glance["reassuring"] += reassuring
    if "pgs" in more:
        notable, reassuring = pgs_mod.bullets(more["pgs"])
        glance["notable"] += notable
        glance["reassuring"] += reassuring
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
    extra = dict(clinvar=cv, clinvar_text=cv_text, pharmcat=pc, pharmcat_text=pc_text, glance=glance, narrative=narrative, pop=pop, card=card,
                 sections=sections)
    rep = write(out / "REPORT.md", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, extra)
    html = write_html(out / "REPORT.html", loaded, qc, ct, roh, roh_text, panel, combo, y, mt, ancestry, input_name=inp.name, extra=extra)
    json.dump(dict(input=str(inp.resolve()), args=[str(Path(x).resolve()) if Path(x).exists() else x for x in run_args],
                   report_html=str((out / "REPORT.html").resolve()), report_md=str((out / "REPORT.md").resolve()), glance=glance,
                   source=loaded.source, build=loaded.build, count=loaded.count, qc=qc,
                   roh_total_mb=float(roh.length_mb.sum()) if not roh.empty else 0.0, y=y, mt=mt,
                   combos=combo, ancestry=ancestry, notes=loaded.notes, kind=loaded.kind,
                   structural=({**more["sv"]["stats"], "findings": more["sv"]["findings"].drop(columns=["rank"]).to_dict("records") if not more["sv"]["findings"].empty else []}
                               if "sv" in more else None),
                   callers={k: {kk: vv for kk, vv in v.items() if kk != "section"} for k, v in called.items()},
                   polygenic=({k: v for k, v in more["pgs"].items() if k != "section"} if "pgs" in more else None),
                   sections=[dict(id=x["id"], title=x["title"], text=x.get("text"), columns=x.get("columns"), rows=x.get("rows"), notes=x.get("notes")) for x in sections],
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
