import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_relationship_estimates():
    from gnomancer.family import estimate_relationship, _phase
    assert estimate_relationship(3500, 0)[0] == "parent / child"
    assert estimate_relationship(2600, 600)[0] == "full siblings"
    assert estimate_relationship(1750, 0)[0].startswith("grandparent")
    assert estimate_relationship(0, 0)[0].startswith("unrelated")
    assert _phase("AG", "AA", "GG") == "A from first parent, G from second"
    assert _phase("AG", "AG", "AG") == "either way"
    assert _phase("GG", "AA", "AG") == "inconsistent with both parents"


def test_mendelian_errors():
    from gnomancer.family import mendelian_errors
    idx = ["rs1", "rs2", "rs3", "rs4"]
    child = pd.DataFrame({"chrom": "1", "pos": [1, 2, 3, 4], "genotype": ["AG", "CC", "TT", None]}, index=idx)
    p1 = pd.DataFrame({"chrom": "1", "pos": [1, 2, 3, 4], "genotype": ["AA", "CT", "AA", "GG"]}, index=idx)
    p2 = pd.DataFrame({"chrom": "1", "pos": [1, 2, 3, 4], "genotype": ["GG", "CC", "TT", "GG"]}, index=idx)
    errs, n = mendelian_errors(child, p1, p2)
    assert (errs, n) == (1, 3)   # rs3: TT impossible with an AA parent; rs4 skipped (no call)


@pytest.mark.slow
def test_family_end_to_end(tmp_path):
    """Synthetic trio: parent/child relationships, sibling detection, no false matches, Mendelian consistency."""
    subprocess.run([sys.executable, str(ROOT / "examples" / "make_family.py"), str(tmp_path / "fam")], check=True, capture_output=True)
    kits = [str(tmp_path / "fam" / f"{n}_AncestryDNA.txt") for n in ("dad", "mom", "child1", "stranger")]
    r = subprocess.run([sys.executable, "-m", "gnomancer", "family", *kits, "--names", "dad,mom,child1,stranger", "--out", str(tmp_path / "out"), "--no-haplogroups"],
                       capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode == 0, r.stderr[-800:]
    md = (tmp_path / "out" / "family-report.md").read_text()
    assert "| dad & child1 |" in md and "**parent / child**" in md
    assert "| dad & mom | 0 |" in md and "| child1 & stranger | 0 |" in md
    assert "Mendelian inconsistencies" in md and "(0.00%)" in md
    assert (tmp_path / "out" / "family-report.html").exists()
