import pandas as pd
import pytest

from gnomancer.load import load
from gnomancer.qc import chromosome_table, summary
from gnomancer.roh import find_roh, interpret
from gnomancer.panel import evaluate, combos, load_markers
from gnomancer.report import write
from gnomancer.report_html import write_html
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize("name,source", [("demo_AncestryDNA.txt", "AncestryDNA"), ("demo_23andMe.txt", "23andMe"),
                                         ("demo_MyHeritage.csv", "MyHeritage")])
def test_load_every_vendor_layout(demo_dir, name, source):
    L = load(str(demo_dir / name))
    assert L.source == source
    assert L.build == 37
    assert set(L.df.chrom) >= {"1", "22", "X", "MT"}
    mt = L.df[L.df.chrom == "MT"].genotype.dropna()
    assert mt.str.len().max() == 1, "MT must be haploid after normalization"


def test_qc_and_sex(demo_dir):
    L = load(str(demo_dir / "demo_AncestryDNA.txt"))
    ct = chromosome_table(L.df)
    qc = summary(L.df, ct, L.sex)
    assert qc["sex_inferred"] == "Female"
    assert 0.2 < qc["autosomal_het"] < 0.4


def test_roh_detects_a_long_homozygous_run():
    pos = list(range(1_000_000, 21_000_000, 10_000))  # 2000 SNPs over 20 Mb
    df = pd.DataFrame({"chrom": "5", "pos": pos, "genotype": ["AA"] * len(pos)}, index=[f"rs{i}" for i in range(len(pos))])
    df.iloc[500, 2] = "AG"  # one tolerated het
    r = find_roh(df)
    assert len(r) == 1 and r.iloc[0].length_mb > 19
    assert "related" in interpret(pd.concat([r] * 12)).lower()


def _df(**gt):
    return pd.DataFrame({"chrom": ["1"] * len(gt), "pos": range(1, len(gt) + 1), "genotype": list(gt.values())}, index=list(gt.keys()))


def test_panel_counts_and_strand_flip():
    markers = [dict(rsid="rs1", category="trait", gene="G", label="x", ref="A", effect="G", interp0="", interp1="one", interp2="", note=""),
               dict(rsid="rs2", category="trait", gene="G", label="y", ref="A", effect="T", interp0="", interp1="", interp2="", note="")]
    out = evaluate(_df(rs1="CT", rs2="AA"), markers)   # rs1 reported on the reverse strand (CT == AG flipped)
    r1 = out[out.rsid == "rs1"].iloc[0]
    assert r1.effect_count == "1/2" and "strand-flipped" in r1.interpretation
    out = evaluate(_df(rs1="GG", rs2="GG"), markers)   # rs2 observed G is impossible for an A/T SNP
    assert "mismatch" in out[out.rsid == "rs2"].iloc[0].interpretation
    assert out[out.rsid == "rs1"].iloc[0].effect_count == "2/2"


def test_combos_apoe_and_abo():
    c = {k: (v, n) for k, v, n in combos(_df(rs429358="TC", rs7412="CC", rs8176719="DD", rs8176746="GG"))}
    assert c["APOE"][0].startswith("e3/e4")
    assert c["ABO blood group"][0].startswith("O")
    c = {k: (v, n) for k, v, n in combos(_df(rs429358="TT", rs7412="CT", rs8176719="II", rs8176746="GT"))}
    assert c["APOE"][0].startswith("e2/e3")
    assert c["ABO blood group"][0].startswith("AB")


def test_reports_render(demo_dir, tmp_path):
    L = load(str(demo_dir / "demo_AncestryDNA.txt"))
    ct = chromosome_table(L.df)
    qc = summary(L.df, ct, L.sex)
    roh = find_roh(L.df)
    panel = evaluate(L.df, load_markers(ROOT / "gnomancer" / "data" / "markers.tsv"))
    y, mt = dict(error="skipped"), dict(error="skipped")
    md = write(tmp_path / "r.md", L, qc, ct, roh, interpret(roh), panel, combos(L.df), y, mt, None)
    html = write_html(tmp_path / "r.html", L, qc, ct, roh, interpret(roh), panel, combos(L.df), y, mt, None, "demo")
    assert "## Marker panel" in md.read_text()
    t = html.read_text()
    assert "<svg" in t and "Marker panel" in t and "<script" in t and "http" not in t.split("<footer>")[0].replace("https://github.com", "")


def test_versions_agree():
    """The version is declared in four files; a release bumps all of them."""
    import json, re
    from gnomancer import __version__
    root = Path(__file__).resolve().parent.parent
    assert re.search(r'^version = "([^"]+)"', (root / "pyproject.toml").read_text(), re.M).group(1) == __version__
    sj = json.loads((root / "server.json").read_text())
    assert sj["version"] == __version__ and sj["packages"][0]["version"] == __version__
    assert json.loads((root / ".claude-plugin" / "plugin.json").read_text())["version"] == __version__
