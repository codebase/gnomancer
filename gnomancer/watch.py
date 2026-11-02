"""`gnomancer watch <report_dir>`: re-check ClinVar for a report and say only what changed.

ClinVar reclassifies variants every month. Clinicians re-check patients' variants by hand; nobody does it for
consumers. This keeps a snapshot of every ClinVar-classified allele the user carries (any significance) next to
the report, downloads a new ClinVar release when NCBI publishes one, rescans the same raw file, and diffs.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import urllib.request
from pathlib import Path

from . import clinvar as cv
from .load import load
from .setup_cmd import CLINVAR_URL, _dl

PATH = ("Pathogenic", "Likely_pathogenic")


def _is_path(sig: str) -> bool:
    return any(s in sig for s in PATH) and not sig.startswith(("Benign", "Likely_benign"))


def snapshot_path(report_dir: Path) -> Path:
    return report_dir / "clinvar-snapshot.json"


def write_snapshot(report_dir: Path, res: dict, input_path: str) -> None:
    json.dump(dict(release=res.get("release"), taken=dt.datetime.now().isoformat(timespec="seconds"), input=input_path,
                   vcf=res.get("vcf"), carried=res.get("carried") or {}), open(snapshot_path(report_dir), "w"))


def remote_last_modified() -> str | None:
    try:
        req = urllib.request.Request(CLINVAR_URL, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as r:  # noqa: S310
            return r.headers.get("Last-Modified")
    except Exception:  # noqa: BLE001
        return None


def diff(old: dict, new: dict) -> dict:
    out = dict(new_pathogenic=[], no_longer_pathogenic=[], review_changed=[], new_records=0, removed_records=0)
    for k, n in new.items():
        o = old.get(k)
        if o is None:
            out["new_records"] += 1
            if _is_path(n["sig"]):
                out["new_pathogenic"].append(dict(key=k, **n, previous=None))
            continue
        if _is_path(n["sig"]) and not _is_path(o["sig"]):
            out["new_pathogenic"].append(dict(key=k, **n, previous=o["sig"]))
        elif _is_path(o["sig"]) and not _is_path(n["sig"]):
            out["no_longer_pathogenic"].append(dict(key=k, **n, previous=o["sig"]))
        elif _is_path(n["sig"]) and o.get("rev") != n.get("rev"):
            out["review_changed"].append(dict(key=k, **n, previous=o.get("rev")))
    out["removed_records"] = sum(1 for k in old if k not in new)
    return out


def render(d: dict, old_release, new_release) -> str:
    L = [f"# ClinVar watch\n\nPrevious release: {old_release or 'unknown'} → current: {new_release or 'unknown'}\n"]
    if not (d["new_pathogenic"] or d["no_longer_pathogenic"] or d["review_changed"]):
        L.append("No change in pathogenic or likely-pathogenic classification at any site you carry.")
    def row(x, extra):
        return f"- {x['gene'] or '?'} {x['rsid']} (genotype {x['gt']}): now **{x['sig'].replace('_', ' ')}** [{x['rev'].replace('_', ' ')}]{extra}"
    if d["new_pathogenic"]:
        L.append("\n## Newly pathogenic / likely pathogenic at a site you carry\n")
        L += [row(x, f"; was {x['previous'].replace('_', ' ')}" if x["previous"] else "; new ClinVar record") for x in d["new_pathogenic"]]
        L.append("\nArray calls at rare sites are often wrong; treat each as 'confirm clinically', not as a diagnosis.")
    if d["no_longer_pathogenic"]:
        L.append("\n## Downgraded (no longer pathogenic / likely pathogenic)\n")
        L += [row(x, f"; was {x['previous'].replace('_', ' ')}") for x in d["no_longer_pathogenic"]]
    if d["review_changed"]:
        L.append("\n## Review status changed on a pathogenic site you carry\n")
        L += [row(x, f"; was [{(x['previous'] or '').replace('_', ' ')}]") for x in d["review_changed"]]
    L.append(f"\n{d['new_records']} new ClinVar record(s) at sites you carry; {d['removed_records']} removed.")
    return "\n".join(L) + "\n"


def run(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="gnomancer watch", description="Re-check ClinVar for a report and report what changed.")
    ap.add_argument("report_dir", help="a gnomancer output folder (contains summary.json)")
    ap.add_argument("--force", action="store_true", help="rescan even if ClinVar has not changed")
    ap.add_argument("--no-download", action="store_true", help="use the ClinVar file already on disk")
    a = ap.parse_args(argv)
    rd = Path(a.report_dir).expanduser()
    summ = rd / "summary.json"
    if not summ.exists():
        print(f"{rd} has no summary.json", file=sys.stderr)
        return 2
    s = json.loads(summ.read_text())
    inp = s.get("input")
    if not inp or not Path(inp).exists():
        print(f"raw file not found at {inp}; cannot rescan", file=sys.stderr)
        return 2
    snap = json.loads(snapshot_path(rd).read_text()) if snapshot_path(rd).exists() else None
    vcf = cv.DEFAULT_VCF
    if not a.no_download:
        lm = remote_last_modified()
        local = dt.datetime.fromtimestamp(vcf.stat().st_mtime, dt.timezone.utc) if vcf.exists() else None
        remote = dt.datetime.strptime(lm, "%a, %d %b %Y %H:%M:%S %Z").replace(tzinfo=dt.timezone.utc) if lm else None
        if remote and (local is None or remote > local):
            print(f"new ClinVar release ({lm}); downloading ...", file=sys.stderr)
            if vcf.exists():
                vcf.unlink()
            _dl(CLINVAR_URL, vcf, "clinvar")
        elif snap and not a.force:
            print(f"ClinVar unchanged since {snap.get('taken')} (release {snap.get('release')}); nothing to do. Use --force to rescan anyway.")
            return 0
    print("rescanning ...", file=sys.stderr)
    df = load(inp).df
    res = cv.scan(df, vcf)
    if snap is None:
        write_snapshot(rd, res, inp)
        print(f"Baseline snapshot written ({len(res['carried']):,} ClinVar-classified alleles carried, release {res.get('release')}). Run again after the next ClinVar release to see changes.")
        return 0
    d = diff(snap["carried"], res["carried"])
    text = render(d, snap.get("release"), res.get("release"))
    (rd / "clinvar-watch.md").write_text(text)
    write_snapshot(rd, res, inp)
    print(text)
    return 0
