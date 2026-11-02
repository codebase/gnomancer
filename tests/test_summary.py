import pandas as pd

from gnomancer.summary import at_a_glance, build_prompt, POP_NAMES


def test_glance_and_prompt_build():
    panel = pd.DataFrame([
        dict(rsid="rs1", category="health", gene="TCF7L2", label="T2D", genotype="TT", effect_count="2/2", interpretation="two risk alleles"),
        dict(rsid="rs2", category="carrier", gene="CFTR", label="dF508", genotype="DI", effect_count="1/2", interpretation="CF carrier"),
        dict(rsid="rs3", category="health", gene="F5", label="Leiden", genotype="CC", effect_count="0/2", interpretation="negative"),
        dict(rsid="rs4", category="pgx", gene="CYP2C9", label="*2", genotype="not on array", effect_count="", interpretation=""),
    ])
    combos = [("APOE", "e3/e4 (rs429358 TC, rs7412 CC)", "one e4 copy"), ("ABO blood group", "O (rs8176719 DD, rs8176746 GG)", "")]
    y = dict(haplogroup="R-M269", ycc="R1b1a2", yleaf=dict(trees=dict(yfull=dict(haplogroup="R-L21*"))))
    mt = dict(haplogroup="H1", quality=1.0)
    anc = dict(population_code="IBS", superpopulation_code="EUR", superpopulation_percent=0.8, population_percent=0.3)
    cv = dict(available=True, hits=pd.DataFrame([dict(gene="ABCA4", zygosity="heterozygous", condition="Stargardt; x", confidence="high"),
                                                   dict(gene="NF1", zygosity="homozygous", condition="NF1", confidence="low")]))
    pc = dict(calls=pd.DataFrame([dict(gene="CYP2C19", diplotype="*1/*2", phenotype="Intermediate Metabolizer"),
                                  dict(gene="TPMT", diplotype="*1/*1", phenotype="Normal Metabolizer")]))
    qc = dict(no_call_rate=0.002, sex_inferred="Male", flags=[])
    g = at_a_glance(None, qc, pd.DataFrame(), "No runs of homozygosity: fine", combos, y, mt, anc, panel, cv, pc)
    text = " ".join(sum(g.values(), []))
    assert "R-M269" in text and "R-L21" in text and POP_NAMES["IBS"] in text
    assert "TCF7L2" in text and "CFTR" in text and "ABCA4" in text and "CYP2C19 *1/*2" in text
    assert "APOE e3/e4" in text and "1 other ClinVar" in text and "TPMT" in text
    prompt = build_prompt(g, combos, y, mt, anc, panel, cv, pc, qc)
    assert "DATA (JSON)" in prompt and "\"rsid\"" not in prompt and "TCF7L2" in prompt
