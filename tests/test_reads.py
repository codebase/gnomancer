"""Read-level analyses: parsing each caller's output, the reference check, mitochondrial merging, polygenic scores.

The callers themselves (ExpansionHunter, PyPGx, T1K, SMNCopyNumberCaller, samtools) are not run here; these tests
feed their output formats to the parsers and the interpretation that follows.
"""
import json

import pytest

from gnomancer import bam, callers, pgs
from gnomancer.seqsites import ASSUMED, READ, Site, genotype, named_positions
from gnomancer.seqvcf import key, read

from test_seq import vcf

# ------------------------------------------------------------------------------------------- repeat expansions

EH = """##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr4\t3074876\t.\tC\t<STR17>,<STR44>\t.\tPASS\tEND=3074933;REF=19;RL=57;RU=CAG;VARID=HTT;REPID=HTT;STR_STATUS=full_mutation;STR_NORMAL_MAX=34;STR_PATHOLOGIC_MIN=36;DisplayRU=CAG;Disease=Huntington;InheritanceMode=AD\tGT:SO:REPCN:REPCI\t1/2:SPANNING/INREPEAT:17/44:17-17/40-52
chr4\t3074939\t.\tG\t.\t.\tPASS\tEND=3074966;REF=9;RL=27;RU=CCG;VARID=HTT_CCG;REPID=HTT_CCG\tGT:SO:REPCN:REPCI\t0/0:SPANNING/SPANNING:9/9:9-9/9-9
chr9\t69037286\t.\tA\t<STR9>,<STR70>\t.\tPASS\tEND=69037304;REF=6;RU=GAA;VARID=FXN;REPID=FXN;STR_STATUS=full_mutation;STR_NORMAL_MAX=35;STR_PATHOLOGIC_MIN=51;DisplayRU=GAA;Disease=FRDA;InheritanceMode=AR\tGT:SO:REPCN:REPCI\t1/2:SPANNING/INREPEAT:9/70:9-9/60-90
chrX\t147912050\t.\tC\t<STR70>\t.\tPASS\tEND=147912110;REF=20;RU=CGG;VARID=FMR1;REPID=FMR1;STR_STATUS=pre_mutation;STR_NORMAL_MAX=55;STR_PATHOLOGIC_MIN=200;DisplayRU=CGG;Disease=FXS;InheritanceMode=XR\tGT:SO:REPCN:REPCI\t1:INREPEAT:70:62-81
chr19\t45770204\t.\tC\t<STR11>,<STR13>\t.\tPASS\tEND=45770264;REF=20;RU=CTG;VARID=DMPK;REPID=DMPK;STR_STATUS=normal;STR_NORMAL_MAX=34;STR_PATHOLOGIC_MIN=50;DisplayRU=CTG;Disease=DM1;InheritanceMode=AD\tGT:SO:REPCN:REPCI\t1/2:SPANNING/SPANNING:11/13:11-11/13-13
chr4\t39348424\t.\tA\t<STR52>\t.\tPASS\tEND=39348483;REF=11;RU=AARRG;VARID=RFC1;REPID=RFC1;STR_STATUS=full_mutation;STR_NORMAL_MAX=11;STR_PATHOLOGIC_MIN=12;DisplayRU=AARRG;Disease=CANVAS;InheritanceMode=AR\tGT:SO:REPCN:REPCI\t0/1:SPANNING/INREPEAT:11/52:11-11/45-74
"""


def test_repeats_are_read_scored_and_ordered_by_concern():
    rows = callers.parse_repeats(EH, "Male")
    by = {r["locus"]: r for r in rows}
    assert "HTT_CCG" not in by                                       # flanking unit: sized, but not what the thresholds mean
    assert [r["locus"] for r in rows] == ["HTT", "FXN", "FMR1", "RFC1", "DMPK"]
    assert by["HTT"]["level"] == "expanded" and by["HTT"]["counts"] == [17, 44] and by["HTT"]["estimated"]
    assert by["HTT"]["disease"] == "Huntington disease"
    assert by["FXN"]["level"] == "carrier" and "recessive" in by["FXN"]["reading"]   # one expanded allele of a recessive locus
    assert by["FMR1"]["level"] == "intermediate" and by["FMR1"]["counts"] == [70]    # one X in a male
    assert by["DMPK"]["level"] == "normal" and not by["DMPK"]["estimated"]
    assert by["RFC1"]["level"] == "uncertain" and "AAGGG" in by["RFC1"]["reading"]   # length alone decides nothing here


def test_repeat_summary_names_only_what_is_actionable():
    rows = callers.parse_repeats(EH, "Male")
    notable, good = callers.repeats_bullets(dict(available=True, loci=rows))
    text = " ".join(notable)
    assert "HTT" in text and "FXN" in text and "FMR1 (70)" in text and "RFC1" not in text and not good
    section = callers.repeats_section(rows)
    assert section["id"] == "repeats" and section["rows"][0][0] == "HTT" and "17/44 (estimated)" in section["rows"][0][3]
    normal = [r for r in rows if r["level"] == "normal"]
    assert callers.repeats_bullets(dict(available=True, loci=normal)) == ([], ["1 of 1 disease-associated repeat loci checked are in the normal range (Huntington, fragile X, myotonic dystrophy and the ataxias among them)."])
    assert callers.repeats_bullets(dict(available=False)) == ([], [])


# ---------------------------------------------------------------------------------------------- CYP2D6, HLA, SMN

def test_pypgx_table_is_parsed_with_or_without_a_sample_header():
    text = "\tGenotype\tPhenotype\tHaplotype1\tHaplotype2\tAlternativePhase\tVariantData\tCNV\nHG003\t*4/*5\tPoor Metabolizer\t*4;\t*5;\t;\t...\tWholeDel1\n"
    r = callers.parse_pypgx(text)
    assert (r["genotype"], r["phenotype"], r["cnv"]) == ("*4/*5", "Poor Metabolizer", "WholeDel1")
    r = callers.parse_pypgx("Genotype\tPhenotype\tHaplotype1\tHaplotype2\tAlternativePhase\tVariantData\tCNV\nS\t*1/*2x2\tUltrarapid Metabolizer\t*1;\t*2;\t;\t...\tWholeDup1\n")
    assert r["genotype"] == "*1/*2x2" and r["cnv"] == "WholeDup1"
    assert callers.parse_pypgx("") is None
    section = callers.cyp2d6_section([dict(gene="CYP2D6", genotype="*4/*5", phenotype="Poor Metabolizer", cnv="WholeDel1")])
    assert section["rows"] == [["CYP2D6", "*4/*5", "Poor Metabolizer", "one copy deleted (*5)"]]
    poor = dict(available=True, calls=[dict(gene="CYP2D6", genotype="*4/*4", phenotype="Poor Metabolizer", cnv="Normal")])
    normal = dict(available=True, calls=[dict(gene="CYP2D6", genotype="*1/*2", phenotype="Normal Metabolizer", cnv="Normal")])
    assert callers.cyp2d6_bullets(poor)[0] and not callers.cyp2d6_bullets(poor)[1]
    assert callers.cyp2d6_bullets(normal)[1] and not callers.cyp2d6_bullets(normal)[0]


T1K = """HLA-A\t2\tHLA-A*01:01:01:01\t120.5\t60\tHLA-A*31:01:02:01\t98.1\t48\t.
HLA-B\t2\tHLA-B*57:01:01:01\t80.0\t55\tHLA-B*27:05:02:01\t77.0\t41\t.
HLA-C\t1\tHLA-C*06:02:01:01\t190.2\t0\t.\t0\t-1\t.
HLA-DQB1\t2\tHLA-DQB1*02:01:01:01\t60.0\t33\tHLA-DQB1*03:02:01:01\t58.0\t30\t.
HLA-DRB1\t2\tHLA-DRB1*03:01:01:01,HLA-DRB1*03:01:01:02\t60.0\t33\tHLA-DRB1*07:01:01:01\t58.0\t30\t.
HLA-E\t2\tHLA-E*01:01:01:01\t60.0\t33\tHLA-E*01:03:02:01\t58.0\t30\t.
"""


def test_hla_types_are_reduced_to_two_fields_and_matched_to_known_alleles():
    types = callers.parse_t1k(T1K)
    assert types["HLA-A"] == [("A*01:01", 60), ("A*31:01", 48)]
    assert types["HLA-C"] == [("C*06:02", 0), ("C*06:02", 0)]                # listed once: homozygous
    assert types["HLA-DRB1"][0][0] == "DRB1*03:01" and "HLA-E" not in types   # ambiguity list collapsed; non-classical genes left out
    notes = {a: (n, m) for a, n, m in callers.hla_notes(types)}
    assert set(notes) == {"A*31:01", "B*57:01", "B*27", "C*06:02", "DRB1*03:01"}
    assert notes["C*06:02"][0] == 2 and "abacavir" in notes["B*57:01"][1]
    notable, good = callers.hla_bullets(dict(available=True, types=types))
    assert any("B*57:01" in n and "abacavir" in n for n in notable) and any("A*31:01" in n for n in notable) and not good
    section = callers.hla_section(types)
    assert ["HLA-C", "C*06:02", "C*06:02", "homozygous", "low confidence"] in section["rows"]
    clean = callers.parse_t1k("HLA-A\t2\tHLA-A*02:01:01:01\t1\t60\tHLA-A*03:01:01:01\t1\t60\t.\n")
    assert callers.hla_bullets(dict(available=True, types=clean)) == ([], ["No HLA allele linked to severe drug reactions (B*57:01, B*15:02, B*58:01, A*31:01) was found."])


def test_hla_prefix_matching_does_not_confuse_neighbouring_alleles():
    types = dict.fromkeys(callers.HLA_GENES, [("X", 1), ("X", 1)])
    types["HLA-B"] = [("B*15:01", 50), ("B*57:03", 50)]                       # neither B*15:02 nor B*57:01
    assert callers.hla_notes(types) == []


@pytest.mark.parametrize("smn1,expected", [(2, "not a carrier"), (1, "carrier of spinal muscular atrophy"), (0, "No functional SMN1"), (3, "3 SMN1 copies")])
def test_smn_copy_number_readings(smn1, expected):
    text = ("Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tSMN2delta7-8_CN\tTotal_CN_raw\tFull_length_CN_raw\tg.27134T>G_CN\n"
            f"S\t{smn1 == 0}\t{smn1 == 1}\t{smn1}\t2\t0\t4.1\t4.0\t0\n")
    r = callers.parse_smn(text)
    assert r["smn1"] == smn1 and r["smn2"] == 2 and r["is_carrier"] is (smn1 == 1)
    assert expected in callers.smn_section(r)["text"]
    notable, good = callers.smn_bullets(dict(available=True, call=r))
    assert bool(notable) == (smn1 < 2) and bool(good) == (smn1 >= 2)


def test_smn_silent_carrier_marker_is_mentioned_and_a_no_call_is_not_reassuring():
    text = "Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\tSMN2delta7-8_CN\tTotal_CN_raw\tFull_length_CN_raw\tg.27134T>G_CN\nS\tFalse\tFalse\t2\t1\t0\t3\t3\t1\n"
    assert "silent" in callers.smn_section(callers.parse_smn(text))["text"]
    nocall = callers.parse_smn("Sample\tisSMA\tisCarrier\tSMN1_CN\tSMN2_CN\nS\tNone\tNone\tNone\tNone\n")
    assert nocall["smn1"] is None and callers.smn_bullets(dict(available=True, call=nocall)) == ([], [])
    assert "could not be called" in callers.smn_section(nocall)["text"]


def test_a_missing_tool_is_reported_not_raised(tmp_path, monkeypatch):
    monkeypatch.setattr(callers.wgstools, "tool_bin", lambda *a, **k: None)
    b = bam.Bam(path=tmp_path / "x.bam", cram=False, build=38, contigs={"1": ("chr1", 248956422)}, reference=None)
    for name, fn, bullets in callers.ALL:
        r = fn(b, tmp_path, "Male") if name == "repeats" else fn(b, tmp_path)
        assert r["available"] is False and "setup --wgs" in r["error"] and bullets(r) == ([], [])


# --------------------------------------------------------------------------------- the reads as a second witness

SITE = Site("rsA", "panel", "1", 5000, "A", ("G",))


def fake_bam(tmp_path):
    return bam.Bam(path=tmp_path / "x.bam", cram=False, build=38, contigs={"1": ("chr1", 248956422), "MT": ("chrM", 16569)}, reference=tmp_path / "ref.fa")


@pytest.mark.parametrize("pile,expected", [
    ((30, ["A"], [30]), ("AA", READ)),
    ((30, ["A", "G"], [29, 1]), ("AA", READ)),                                   # a stray read is not a variant
    ((3, ["A"], [3]), (None, "covered too thinly to call")),
    ((30, ["A", "G"], [15, 15]), (None, "the reads show a variant here that the VCF does not list")),
    (None, (None, "covered too thinly to call")),                                # no reads at all
])
def test_silent_positions_are_settled_by_the_reads(tmp_path, monkeypatch, pile, expected):
    p = vcf(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")])
    named = named_positions([SITE])
    d = read(p, {key("1", 5000)}, named)
    assert genotype(d, SITE, True) == ("AA", ASSUMED)
    monkeypatch.setattr(bam, "pileup_sites", lambda b, sites: {("1", 5000): pile} if pile else {("1", 1): (30, ["A"], [30])})
    stats = bam.confirm_reference(fake_bam(tmp_path), d, named)
    assert d.reads_checked and d.observed and stats["checked"] == 1
    assert genotype(d, SITE, True) == expected
    assert stats["confirmed"] + stats["thin"] + stats["disagree"] == 1


def test_positions_with_a_call_are_not_rechecked(tmp_path, monkeypatch):
    p = vcf(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")])
    named = named_positions([SITE])
    d = read(p, {key("1", 5000)}, named)
    asked = []
    monkeypatch.setattr(bam, "pileup_sites", lambda b, sites: asked.extend(sites) or {})
    assert bam.confirm_reference(fake_bam(tmp_path), d, named)["checked"] == 0 and asked == []


def test_reads_can_flag_a_sample_mismatch(tmp_path):
    p = vcf(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")])
    d = read(p, set(), {})
    d.sex = "Female"
    b = fake_bam(tmp_path)
    b.cov = dict(mean_depth=9.0, sex="Male")
    b.confirm = dict(checked=1000, disagree=80, examples=["chr1:5", "chr2:9"])
    flags = " ".join(bam.flags(b, d))
    assert "9x" in flags and "may not be from the same person" in flags and "may not be from the same sample" in flags
    b.cov, b.confirm = dict(mean_depth=34.0, sex="Female"), dict(checked=1000, disagree=1, examples=[])
    assert bam.flags(b, d) == []


def test_read_level_mitochondrial_variants_join_the_call_set_once(tmp_path):
    p = vcf(tmp_path, [("chrM", 263, "A", "G", "1/1:0,2000:2000:99")])
    d = read(p, set(), {})
    het = dict(variants=[dict(pos=263, ref="A", alt="G", share=0.999, depth=2000, forward=900, reverse=1100, homoplasmic=True),
                         dict(pos=3243, ref="A", alt="G", share=0.04, depth=9000, forward=200, reverse=160, homoplasmic=False)])
    assert bam.merge_mito(d, het) == 1
    (c,) = d.calls[key("MT", 3243)]
    assert c.ok and c.ab == 0.04 and len(d.calls[key("MT", 263)]) == 1
    assert bam.merge_mito(d, het) == 0


def test_mitochondrial_artifact_positions_are_known():
    assert {310, 3107, 16183} <= bam.MT_ARTIFACT and 3243 not in bam.MT_ARTIFACT and 263 not in bam.MT_ARTIFACT


# ------------------------------------------------------------------------------------------ polygenic scores

@pytest.fixture
def tiny_pgs(tmp_path, monkeypatch):
    d = tmp_path / "pgs"
    d.mkdir()
    (d / "PGS1.tsv").write_text("# test\nchrom\tpos\tref\talt\teffect\tweight\taf_AFR\taf_AMR\taf_EAS\taf_EUR\taf_SAS\n"
                                "1\t5000\tA\tG\tG\t1.0\t0.5\t0.5\t0.5\t0.5\t0.5\n"       # effect allele is the alternate
                                "1\t6000\tC\tT\tC\t0.5\t0.5\t0.5\t0.5\t0.5\t0.5\n"       # effect allele is the reference
                                "1\t7000\tG\tA\tA\t2.0\t0.25\t0.25\t0.25\t0.25\t0.25\n"
                                "1\t8000\tGT\tG\tG\t3.0\t0.25\t0.25\t0.25\t0.25\t0.25\n")      # a deletion is the effect allele
    ref = {p: dict(n=5, mean=2.0, sd=1.0, scores=[0.0, 1.0, 2.0, 3.0, 4.0]) for p in pgs.POPS}
    (d / "index.json").write_text(json.dumps({"PGS1": dict(id="PGS1", label="Test condition", sex=None, note="", citation="Someone et al., Journal 2020",
                                                           variants_published=3, variants_used=3, reference=ref)}))
    monkeypatch.setattr(pgs, "PGS_DIR", d)
    return d


def score(tmp_path, rows, pop="EUR", sex="Male", whole=True):
    named = pgs.positions()
    d = read(vcf(tmp_path, rows), {key(c, p) for c, ps in named.items() for p in ps}, named)
    return pgs.score(d, whole, sex, pop)


def test_polygenic_score_counts_effect_alleles_including_reference_ones(tmp_path, tiny_pgs):
    # 5000: het for the effect allele (1 x 1.0). 6000: absent, reference C is the effect allele (2 x 0.5). 7000: hom alt (2 x 2.0).
    res = score(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99"), ("chr1", 7000, "G", "A", "1/1:0,30:30:99")])
    (r,) = res["results"]
    assert r["score"] == 6.0 and r["missing"] == 0 and r["percentile"] == 100.0 and r["z"] == 4.0
    with_del = score(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99"), ("chr1", 7000, "G", "A", "1/1:0,30:30:99"), ("chr1", 8000, "GT", "G", "0/1:15,15:30:99")])
    assert with_del["results"][0]["score"] == 9.0                               # one copy of the deletion adds 3.0
    assert res["section"]["rows"][0][:2] == ["Test condition", "99th percentile"] and not res["population_assumed"]
    assert pgs.bullets(res)[0] and not pgs.bullets(res)[1]
    low = score(tmp_path, [("chr1", 6000, "C", "T", "1/1:0,30:30:99")])["results"][0]
    assert low["score"] == 0.0 and low["percentile"] == 10.0                      # tied with the lowest reference sample: mid-rank


def test_polygenic_score_imputes_a_few_missing_sites_and_refuses_too_many(tmp_path, tiny_pgs):
    # The lightest site (weight 0.5 of 6.5 in total, under the 10% limit) fails quality: imputed at 2 x 0.5 x 0.5 and counted.
    res = score(tmp_path, [("chr1", 6000, "C", "T", "0/1:2,2:4:99")])["results"][0]
    assert res["missing"] == 1 and res["usable"] and res["score"] == 0.5
    # The heaviest site failing is too much of the score to guess.
    res = score(tmp_path, [("chr1", 8000, "GT", "G", "0/1:2,2:4:99")])["results"][0]
    assert res["missing"] == 1 and not res["usable"] and "too many" in res["reason"]
    assert score(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")], whole=False)["available"] is False


def test_polygenic_score_respects_sex_and_says_when_the_population_was_assumed(tmp_path, tiny_pgs):
    idx = json.loads((tiny_pgs / "index.json").read_text())
    idx["PGS1"]["sex"] = "Female"
    (tiny_pgs / "index.json").write_text(json.dumps(idx))
    assert score(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")], sex="Male")["results"] == []
    res = score(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")], sex="Female", pop=None)
    assert len(res["results"]) == 1 and res["population"] == "EUR" and res["population_assumed"] and "run with --ancestry" in res["section"]["text"]


def test_percentiles_are_written_as_proper_ordinals():
    assert [pgs.ordinal(p) for p in (0.2, 1, 2, 3, 4, 11, 12, 13, 21.6, 32.4, 53, 99.8, 100)] == \
        ["1st", "1st", "2nd", "3rd", "4th", "11th", "12th", "13th", "22nd", "32nd", "53rd", "99th", "99th"]


def test_percentile_bands():
    assert [pgs.band(p) for p in (99, 85, 50, 10, 2)] == ["in the top 5%", "higher than about four in five people", "in the middle range",
                                                         "lower than about four in five people", "in the bottom 5%"]


def test_celiac_reading_depends_on_the_dq_types_present():
    dq = lambda a1, a2, b1, b2: {"HLA-DQA1": [(a1, 30), (a2, 30)], "HLA-DQB1": [(b1, 30), (b2, 30)]}  # noqa: E731
    assert "very unlikely" in callers.celiac(dq("DQA1*01:02", "DQA1*01:03", "DQB1*06:02", "DQB1*06:03"))
    assert "DQ2.5 present" in callers.celiac(dq("DQA1*05:01", "DQA1*01:02", "DQB1*02:01", "DQB1*06:02"))
    both = callers.celiac(dq("DQA1*02:01", "DQA1*03:01", "DQB1*02:02", "DQB1*03:02"))
    assert "DQ8 and DQ2.2 present" in both and "does not do here" in both
    assert callers.celiac({}) == ""
