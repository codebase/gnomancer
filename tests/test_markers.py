import collections
from pathlib import Path

from gnomancer.panel import load_markers, COMP

ROOT = Path(__file__).resolve().parent.parent


def test_markers_table_is_well_formed():
    rows = load_markers(ROOT / "gnomancer" / "data" / "markers.tsv")
    assert len(rows) > 50
    dup = [k for k, v in collections.Counter(r["rsid"] for r in rows).items() if v > 1]
    assert not dup, dup
    for r in rows:
        assert r["rsid"].startswith("rs"), r
        assert r["ref"] in COMP and r["effect"] in COMP and r["ref"] != r["effect"], r["rsid"]
        assert r["category"] in {"ancestry", "trait", "health", "carrier", "pgx"}, r["rsid"]
        assert r["label"] and r["gene"], r["rsid"]


def test_haplogroup_context_lookup():
    from gnomancer.context import lookup
    names = [r["name"] for r in lookup("Y", "R-M269", "R1b1a2a1a2")]
    assert names[0] == "R-M269" and names[-1] == "Haplogroup R"
    assert lookup("MT", "H1")[0]["name"] == "Haplogroup H1"
    assert lookup("MT", "999") == []


def test_frequencies_table_and_genotype_frequency():
    from gnomancer.panel import load_frequencies, genotype_frequency, evaluate, load_markers
    import pandas as pd
    fr = load_frequencies()
    assert len(fr) >= 90 and abs(fr["rs1426654"]["EUR"] - 0.997) < 0.01
    assert abs(genotype_frequency(2, 2, 0.5) - 0.25) < 1e-9 and abs(genotype_frequency(1, 2, 0.5) - 0.5) < 1e-9
    assert genotype_frequency(1, 1, 0.3) == 0.3 and abs(genotype_frequency(0, 1, 0.3) - 0.7) < 1e-9
    df = pd.DataFrame({"chrom": ["15"], "pos": [1], "genotype": ["AA"]}, index=["rs1426654"])
    out = evaluate(df, [m for m in load_markers(ROOT / "gnomancer" / "data" / "markers.tsv") if m["rsid"] == "rs1426654"])
    assert out.iloc[0]["gf_EUR"] > 0.99 and out.iloc[0]["gf_AFR"] < 0.01
