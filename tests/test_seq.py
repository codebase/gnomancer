"""Sequencing VCFs: reading, quality gates, genotypes at named sites, ClinVar, PharmCAT input, structural variants.

Every fixture here is a few hand-written VCF lines, so each test says exactly what goes in and what must come out.
"""
import gzip

import pandas as pd
import pytest

from gnomancer import clinvar, pharmcat, seqqc, sv
from gnomancer.panel import evaluate
from gnomancer.seqsites import ASSUMED, CALLED, READ, Site, genotype, load_sites, named_positions, view
from gnomancer.seqvcf import fetch, key, norm_chrom, read, sniff, trim

HEADER38 = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n##contig=<ID=chrM,length=16569>\n"
HEADER37 = "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def vcf(tmp_path, rows, header=HEADER38, name="s.vcf", gz=False):
    """rows: (chrom, pos, ref, alt, sample) with optional filter, format, info."""
    lines = []
    for r in rows:
        chrom, pos, ref, alt, sample = r[:5]
        filt = r[5] if len(r) > 5 else "PASS"
        fmt = r[6] if len(r) > 6 else "GT:AD:DP:GQ"
        info = r[7] if len(r) > 7 else "."
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t{info}\t{fmt}\t{sample}\n")
    p = tmp_path / (name + (".gz" if gz else ""))
    text = header + COLS + "".join(lines)
    if gz:
        with gzip.open(p, "wt") as f:
            f.write(text)
    else:
        p.write_text(text)
    return p


def load(tmp_path, rows, sites=(), header=HEADER38, want=(), **kw):
    p = vcf(tmp_path, rows, header)
    named = named_positions(sites)
    w = {key(c, q) for c, ps in named.items() for q in ps} | {key(c, q) for c, q in want}
    return read(p, w, named, **kw)


# ---------------------------------------------------------------------------------------------------- basics

def test_trim_drops_shared_bases_and_keeps_the_left_anchor():
    assert trim(100, "ATCTT", "AT") == (100, "ATCT", "A")           # shared trailing T
    assert trim(100, "CAT", "CATAT") == (100, "C", "CAT")            # insertion written against a longer reference
    assert trim(100, "GA", "GT") == (101, "A", "T")                  # an SNV inside a two-base record moves right
    assert trim(100, "A", "G") == (100, "A", "G")


def test_chromosome_names_are_normalized_and_off_assembly_contigs_dropped():
    assert [norm_chrom(c) for c in ("chr7", "7", "chrM", "MT", "chrX", "chr1_KI270706v1_random", "HLA-A*01:01", "chrEBV")] == \
        ["7", "7", "MT", "MT", "X", None, None, None]


def test_sniff_tells_sequencing_from_an_array_export(tmp_path):
    seq = sniff(vcf(tmp_path, [("chr1", 1000, "A", "G", "0/1:10,12:22:99"), ("chr1", 2000, "AT", "A", "1/1:0,30:30:90")]))
    assert seq.is_sequencing and seq.build == 38 and seq.mt_length == 16569 and not seq.gvcf
    array = tmp_path / "array.vcf"
    array.write_text(HEADER37 + COLS + "".join(f"1\t{1000 + i}\trs{i}\tA\tG\t.\t.\t.\tGT\t0/1\n" for i in range(50)))
    a = sniff(array)
    assert not a.is_sequencing and a.build == 37 and a.rsid_fill == 1.0
    assert not sniff(tmp_path / "missing.vcf").is_sequencing
    not_vcf = tmp_path / "raw.txt"
    not_vcf.write_text("# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAG\n")
    assert not sniff(not_vcf).is_sequencing


def test_sniff_recognises_a_gvcf_and_gzip(tmp_path):
    p = vcf(tmp_path, [("chr1", 1, "N", "<NON_REF>", "0/0:30:90:28", ".", "GT:DP:GQ:MIN_DP", "END=999"),
                       ("chr1", 1000, "A", "G,<NON_REF>", "0/1:10,12,0:22:99")], gz=True)
    assert sniff(p).gvcf and sniff(p).is_sequencing


# ---------------------------------------------------------------------------------------------- quality gates

@pytest.mark.parametrize("sample,filt,ok,why", [
    ("0/1:15,15:30:99", "PASS", True, ""),
    ("0/1:15,15:30:99", ".", True, ""),
    ("0/1:15,15:30:99", "LowQual", False, "caller filter"),
    ("0/1:2,2:4:99", "PASS", False, "only 4 reads"),
    ("0/1:15,15:30:12", "PASS", False, "genotype quality 12"),
    ("0/1:27,3:30:99", "PASS", False, "called heterozygous"),
    ("1/1:10,20:30:99", "PASS", False, "called homozygous"),
    ("1/1:0,30:30:99", "PASS", True, ""),
])
def test_quality_gates(tmp_path, sample, filt, ok, why):
    d = load(tmp_path, [("chr1", 5000, "A", "G", sample, filt)], want=[("1", 5000)])
    (c,) = d.calls[key("1", 5000)]
    assert c.ok is ok and why in c.why
    assert d.stats["passing_alleles"] == int(ok) and d.stats["failed_alleles"] == int(not ok)


def test_a_call_without_evidence_fields_passes_on_the_filter_alone(tmp_path):
    d = load(tmp_path, [("chr1", 5000, "A", "G", "0/1", "PASS", "GT")], want=[("1", 5000)])
    assert d.calls[key("1", 5000)][0].ok


def test_multiallelic_records_are_split_and_each_allele_trimmed(tmp_path):
    d = load(tmp_path, [("chr1", 100, "CAT", "C,CATAT", "1/2:0,14,15:29:99")], want=[("1", 100)])
    calls = {(c.ref, c.alt): c for c in d.calls[key("1", 100)]}
    assert set(calls) == {("CAT", "C"), ("C", "CAT")}
    assert all(c.copies == 1 and c.other_alt and c.zygosity == "heterozygous" for c in calls.values())


def test_mitochondrial_low_share_is_heteroplasmy_not_a_failed_call(tmp_path):
    d = load(tmp_path, [("chrM", 3243, "A", "G", "0/1:1600,400:2000:99")])
    (c,) = d.calls[key("MT", 3243)]
    assert c.ok and c.ab == pytest.approx(0.2)


def test_y_and_mt_calls_are_all_kept_and_other_calls_only_when_asked_for(tmp_path):
    d = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99"), ("chrY", 5000000, "C", "T", "1:0,15:15:99"),
                        ("chrM", 263, "A", "G", "1/1:0,900:900:99"), ("chrUn_KI270742v1", 50, "A", "G", "0/1:15,15:30:99")],
             header=HEADER38 + "##contig=<ID=chrY,length=57227415>\n")
    assert key("1", 100) not in d.calls and key("Y", 5000000) in d.calls and key("MT", 263) in d.calls
    assert d.calls[key("Y", 5000000)][0].zygosity == "hemizygous"
    assert d.stats["off_primary_records"] == 1 and d.stats["snvs"] == 3


def test_build_must_be_known(tmp_path):
    p = vcf(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")], header="##fileformat=VCFv4.2\n")
    with pytest.raises(SystemExit, match="--build"):
        read(p, set(), {})
    assert read(p, set(), {}, build=38).build == 38


# ----------------------------------------------------------------------------------- genotypes at named sites

SNV = Site("rsA", "panel", "1", 5000, "A", ("G",))
DEL = Site("rsDel", "panel", "1", 7000, "GT", ("G",))          # the reference has the longer allele
INS = Site("rsIns", "panel", "1", 9000, "T", ("TC",))          # the reference has the shorter allele (ABO-style)
XSITE = Site("rsX", "panel", "X", 50000000, "C", ("T",))


def gt(d, site, whole=True):
    return genotype(d, site, whole)


def test_snv_site_called_from_a_record(tmp_path):
    d = load(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")], [SNV])
    assert gt(d, SNV) == ("AG", CALLED)
    d = load(tmp_path, [("chr1", 5000, "A", "G", "1/1:0,30:30:99")], [SNV])
    assert gt(d, SNV) == ("GG", CALLED)
    d = load(tmp_path, [("chr1", 5000, "A", "T", "0/1:15,15:30:99")], [SNV])
    assert gt(d, SNV) == ("AT", CALLED)                          # an allele the table does not list is still what was seen


def test_indel_sites_use_array_letters_with_i_for_the_longer_allele(tmp_path):
    d = load(tmp_path, [("chr1", 7000, "GT", "G", "0/1:15,15:30:99"), ("chr1", 9000, "T", "TC", "1/1:0,30:30:99")], [DEL, INS])
    assert gt(d, DEL) == ("ID", CALLED)          # one copy deleted
    assert gt(d, INS) == ("II", CALLED)          # insertion on both copies
    empty = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")], [DEL, INS])
    assert gt(empty, DEL) == ("II", ASSUMED)     # reference = the longer allele = I
    assert gt(empty, INS) == ("DD", ASSUMED)     # reference = the shorter allele = D


def test_indel_written_with_a_shared_trailing_base_still_matches(tmp_path):
    d = load(tmp_path, [("chr1", 7000, "GTA", "GA", "0/1:15,15:30:99")], [DEL])
    assert gt(d, DEL) == ("ID", CALLED)


def test_absent_site_is_reference_only_in_a_whole_genome(tmp_path):
    d = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")], [SNV])
    assert gt(d, SNV, whole=True) == ("AA", ASSUMED)
    g, basis = gt(d, SNV, whole=False)
    assert g is None and "partial" in basis


def test_low_quality_and_missing_calls_are_no_calls_with_a_reason(tmp_path):
    d = load(tmp_path, [("chr1", 5000, "A", "G", "0/1:2,2:4:99")], [SNV])
    g, basis = gt(d, SNV)
    assert g is None and "only 4 reads" in basis
    d = load(tmp_path, [("chr1", 5000, "A", "G", "./.:0,0:0:0")], [SNV])
    assert gt(d, SNV) == (None, "no genotype called here")


def test_a_deletion_spanning_a_site_blocks_the_reference_inference(tmp_path):
    d = load(tmp_path, [("chr1", 4998, "CTAG", "C", "0/1:15,15:30:99")], [SNV])
    g, basis = gt(d, SNV)
    assert g is None and "overlaps" in basis


def test_reference_base_mismatch_is_reported_not_interpreted(tmp_path):
    d = load(tmp_path, [("chr1", 5000, "C", "G", "0/1:15,15:30:99")], [SNV])
    g, basis = gt(d, SNV)
    assert g is None and "wrong build" in basis


def test_gvcf_reference_blocks_decide_between_read_thin_and_not_covered(tmp_path):
    block = lambda start, end, dp, gq: ("chr1", start, "N", "<NON_REF>", f"0/0:{dp}:{gq}:{dp}", ".", "GT:DP:GQ:MIN_DP", f"END={end}")  # noqa: E731
    thin = Site("rsThin", "panel", "1", 12000, "C", ("T",))
    gap = Site("rsGap", "panel", "1", 20000, "C", ("T",))
    d = load(tmp_path, [block(4000, 6000, 30, 90), block(11000, 13000, 3, 9), ("chr1", 30000, "A", "G,<NON_REF>", "0/1:15,15,0:30:99")], [SNV, thin, gap])
    assert d.gvcf
    assert gt(d, SNV) == ("AA", READ)
    assert gt(d, thin) == (None, "covered too thinly to call")
    assert gt(d, gap) == (None, "not covered")
    assert d.stats["reference_blocks"] == 2 and d.stats["snvs"] == 1


def test_male_x_is_one_copy_however_the_caller_writes_it(tmp_path):
    rows = [("chrX", 3000000 + i * 1000, "A", "G", "1/1:0,30:30:99") for i in range(600)]
    d = load(tmp_path, rows + [("chrX", 50000000, "C", "T", "1/1:0,30:30:99")], [XSITE],
             header=HEADER38 + "##contig=<ID=chrX,length=156040895>\n")
    assert d.sex == "Male"
    assert gt(d, XSITE) == ("T", CALLED)
    het = [("chrX", 3000000 + i * 1000, "A", "G", "0/1:15,15:30:99" if i % 2 else "1/1:0,30:30:99") for i in range(600)]
    f = load(tmp_path, het, [XSITE], header=HEADER38 + "##contig=<ID=chrX,length=156040895>\n")
    assert f.sex == "Female" and gt(f, XSITE) == ("CC", ASSUMED)


def test_view_has_the_shape_the_panel_expects(tmp_path):
    d = load(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")], [SNV, DEL])
    df = view(d, {"rsA": SNV, "rsDel": DEL}, True)
    assert list(df.columns) == ["chrom", "pos", "genotype", "basis"]
    assert df.loc["rsA", "genotype"] == "AG" and df.loc["rsDel", "genotype"] == "II"


def test_bundled_sites_cover_the_panel_on_both_builds():
    from gnomancer.panel import load_markers
    from gnomancer.paths import DATA
    markers = {m["rsid"]: m for m in load_markers(DATA / "markers.tsv")}
    for build in (37, 38):
        sites = load_sites(build)
        assert set(markers) <= set(sites)
        for rsid, m in markers.items():
            assert sites[rsid].indel == (m["ref"] in "ID" or m["effect"] in "ID"), rsid
    # BRCA1 185delAG: the deletion, not the duplication ClinVar files under the same rsID.
    assert (load_sites(38)["rs80357914"].ref, load_sites(38)["rs80357914"].alts) == ("ACT", ("A",))
    # ABO: the reference genome is the O (deleted) allele, so the functional allele is written as an insertion.
    assert (load_sites(38)["rs8176719"].ref, load_sites(38)["rs8176719"].alts) == ("T", ("TC",))


def test_grch37_reference_carries_factor_v_leiden_and_absence_is_not_reported_as_a_finding(tmp_path):
    """On GRCh37 the reference base at rs6025 is the Leiden allele. A genome without a record there must not be
    called homozygous Leiden; a normal genome shows a homozygous variant to C, which reads as negative."""
    from gnomancer.panel import load_markers
    from gnomancer.paths import DATA
    f5 = load_sites(37)["rs6025"]
    markers = [m for m in load_markers(DATA / "markers.tsv") if m["rsid"] == "rs6025"]
    assert f5.ref == markers[0]["effect"] == "T"
    absent = load(tmp_path, [("1", 100, "A", "G", "0/1:15,15:30:99")], [f5], header=HEADER37)
    row = evaluate(view(absent, {"rs6025": f5}, True), markers).iloc[0]
    assert row.genotype == "no call" and "reference genome carries" in row.basis
    normal = load(tmp_path, [("1", f5.pos, "T", "C", "1/1:0,30:30:99")], [f5], header=HEADER37)
    row = evaluate(view(normal, {"rs6025": f5}, True), markers).iloc[0]
    assert row.genotype == "CC" and row.effect_count == "0/2"


# -------------------------------------------------------------------------------------------------- ClinVar

CLINVAR = """##fileformat=VCFv4.1
##fileDate=2026-01-01
1\t5000\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=BRCA2:675;CLNDN=Hereditary_breast_ovarian_cancer_syndrome;MC=SO:0001583|missense_variant;RS=111
1\t5000\t2\tA\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=BRCA2:675;CLNDN=Something_else
1\t7000\t3\tGT\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;GENEINFO=CFTR:1080;CLNDN=Inborn_genetic_diseases|Cystic_fibrosis;MC=SO:0001589|frameshift_variant
1\t8000\t4\tC\tT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=ATP7B:540;CLNDN=Wilson_disease
1\t8500\t5\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=ATP7B:540;CLNDN=Wilson_disease
1\t9000\t6\tC\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=PAH:5053;CLNDN=Phenylketonuria
1\t9500\t7\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=XYZ:1
1\t9600\t8\tC\tT\t.\t.\tCLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=XYZ:1
1\t9700\t9\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=HFE:3077;CLNDN=Hemochromatosis_type_1;AF_EXAC=0.08
1\t9800\t10\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided;GENEINFO=GENE9:9;CLNDN=Thing
1\t9900\t11\tG\tA\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNREVSTAT=criteria_provided,_conflicting_classifications;GENEINFO=GENE10:10
MT\t3243\t12\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=MT-TL1:4567;CLNDN=MELAS_syndrome
"""


@pytest.fixture
def clinvar_vcf(tmp_path):
    p = tmp_path / "clinvar_GRCh38.vcf.gz"
    with gzip.open(p, "wt") as f:
        f.write(CLINVAR)
    return p


def scan(tmp_path, clinvar_vcf, rows):
    want = clinvar.positions(38, clinvar_vcf)
    d = read(vcf(tmp_path, rows), want, {})
    return clinvar.scan_seq(d, clinvar_vcf)


def test_clinvar_positions_are_cached_and_cover_every_record(clinvar_vcf):
    pos = clinvar.positions(38, clinvar_vcf)
    assert key("1", 7000) in pos and key("MT", 3243) in pos and len(pos) == 11
    assert clinvar_vcf.with_suffix(".positions.npy").exists()
    assert clinvar.positions(38, clinvar_vcf) == pos


def test_clinvar_matches_exact_alleles_including_indels(tmp_path, clinvar_vcf):
    res = scan(tmp_path, clinvar_vcf, [("chr1", 5000, "A", "G", "0/1:15,15:30:99"),      # pathogenic SNV in an actionable gene
                                       ("chr1", 7000, "GTA", "GA", "0/1:15,15:30:99"),   # pathogenic deletion, untrimmed
                                       ("chr1", 9500, "C", "T", "0/1:15,15:30:99"),      # benign
                                       ("chr1", 9600, "C", "T", "0/1:15,15:30:99")])     # uncertain
    h = res["hits"].set_index("gene")
    assert set(h.index) == {"BRCA2", "CFTR"}
    assert h.loc["BRCA2", "acmg"] and "dominant" in h.loc["BRCA2", "role"] and h.loc["BRCA2", "genotype"] == "A>G"
    assert h.loc["CFTR", "genotype"] == "delT" and h.loc["CFTR", "lof"] and not h.loc["CFTR", "acmg"]
    assert h.loc["CFTR", "condition"].startswith("Cystic fibrosis")          # the catch-all label goes last
    assert res["release"] == "2026-01-01" and res["sites_in_clinvar"] == 4
    assert len(res["carried"]) == 3                                           # benign records are not tracked for the watch


def test_clinvar_same_position_other_allele_is_not_a_match(tmp_path, clinvar_vcf):
    res = scan(tmp_path, clinvar_vcf, [("chr1", 5000, "A", "C", "0/1:15,15:30:99")])
    assert res["hits"].empty and res["sites_in_clinvar"] == 0
    assert "None is classified" in clinvar.interpret_seq(res)


def test_clinvar_groups_weak_calls_common_alleles_and_unreviewed_records(tmp_path, clinvar_vcf):
    res = scan(tmp_path, clinvar_vcf, [("chr1", 9000, "C", "T", "0/1:27,3:30:99"),       # 10% of reads
                                       ("chr1", 9700, "G", "A", "1/1:0,30:30:99"),       # common allele, actionable recessive gene
                                       ("chr1", 9800, "G", "A", "0/1:15,15:30:99"),      # nobody reviewed it
                                       ("chr1", 9900, "G", "A", "0/1:15,15:30:99")])     # conflicting: not a pathogenic record
    h = res["hits"].set_index("gene")
    assert dict(h.confidence) == {"PAH": "low", "HFE": "medium", "GENE9": "low"}
    assert "10% of reads" in h.loc["PAH", "flags"] and res["failed_quality"] == 1
    groups = {k: list(df.gene) for k, _, _, df in clinvar.seq_groups(res["hits"])}
    assert groups == {"actionable": ["HFE"], "aside": ["PAH", "GENE9"]}


def test_clinvar_two_variants_in_one_gene_report_their_phase(tmp_path, clinvar_vcf):
    fmt = "GT:AD:DP:GQ:PS"
    trans = scan(tmp_path, clinvar_vcf, [("chr1", 8000, "C", "T", "0|1:15,15:30:99:7000", "PASS", fmt), ("chr1", 8500, "G", "A", "1|0:15,15:30:99:7000", "PASS", fmt)])
    assert trans["hits"]["flags"].str.contains("on opposite copies").all()
    cis = scan(tmp_path, clinvar_vcf, [("chr1", 8000, "C", "T", "0|1:15,15:30:99:7000", "PASS", fmt), ("chr1", 8500, "G", "A", "0|1:15,15:30:99:7000", "PASS", fmt)])
    assert cis["hits"]["flags"].str.contains("on the same copy").all()
    unphased = scan(tmp_path, clinvar_vcf, [("chr1", 8000, "C", "T", "0/1:15,15:30:99"), ("chr1", 8500, "G", "A", "0/1:15,15:30:99")])
    assert unphased["hits"]["flags"].str.contains("not phased").all()
    assert unphased["hits"].role.str.contains("carrier").all()


def test_clinvar_mitochondrial_hits_are_described_by_read_share(tmp_path, clinvar_vcf):
    res = scan(tmp_path, clinvar_vcf, [("chrM", 3243, "A", "G", "0/1:1600,400:2000:99")])
    assert res["hits"].iloc[0].zygosity == "heteroplasmic (20% of reads)"
    assert [k for k, *_ in clinvar.seq_groups(res["hits"])] == ["both"]
    res = scan(tmp_path, clinvar_vcf, [("chrM", 3243, "A", "G", "1/1:10,1990:2000:99")])
    assert res["hits"].iloc[0].zygosity == "homoplasmic"
    low = scan(tmp_path, clinvar_vcf, [("chrM", 3243, "A", "G", "0/1:9700,300:10000:99")])
    assert [k for k, *_ in clinvar.seq_groups(low["hits"])] == ["mito_low"]        # 3% of reads is not 'on the only copy'


def test_acmg_list_is_v3_3():
    genes = clinvar.acmg_genes()
    assert len(genes) == 84 and {"ABCD1", "CYP27A1", "PLN"} <= set(genes)
    assert {g for g, r in genes.items() if r["inheritance"] == "AR"} == {"MUTYH", "CASQ2", "TRDN", "BTD", "CYP27A1", "GAA", "HFE", "RPE65", "ATP7B"}


# ------------------------------------------------------------------------------------------------- PharmCAT

POSITIONS = """##fileformat=VCFv4.2
##source=test
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t5000\trs1\tA\tG\t.\tPASS\tPX=GENE1
chr1\t5100\trs2\tC\tT\t.\tPASS\tPX=GENE1
chr1\t6000\trs3\tCAT\tC,CATAT\t.\tPASS\tPX=GENE2
chr1\t6500\trs4\tG\tA\t.\tPASS\tPX=GENE2
chr1\t6600\trs5\tG\tA\t.\tPASS\tPX=GENE2
"""


@pytest.fixture
def pc_positions(tmp_path, monkeypatch):
    p = tmp_path / "pharmcat_positions.vcf"
    p.write_text(POSITIONS)
    monkeypatch.setattr(pharmcat, "POSITIONS", p)
    return p


def pc_rows(tmp_path, rows, whole=True, header=HEADER38):
    named = pharmcat.named_keys(pharmcat.positions())
    d = read(vcf(tmp_path, rows, header), {key(c, p) for c, ps in named.items() for p in ps}, named)
    out = tmp_path / "pc.vcf"
    stats = pharmcat.build_vcf_seq(d, out, whole)
    body = {l.split("\t")[2]: l.rstrip("\n").split("\t") for l in out.read_text().splitlines() if not l.startswith("#")}
    return stats, body


def test_pharmcat_input_keeps_indels_phase_and_fills_reference(tmp_path, pc_positions):
    fmt = "GT:AD:DP:GQ:PS"
    stats, body = pc_rows(tmp_path, [("chr1", 5000, "A", "G", "0|1:15,15:30:99:4000", "PASS", fmt), ("chr1", 5100, "C", "T", "1|0:15,15:30:99:4000", "PASS", fmt),
                                     ("chr1", 6000, "C", "CAT", "1/1:0,30:30:99"), ("chr1", 6500, "G", "A", "0/1:2,2:4:99")])
    assert body["rs1"][8:] == ["GT:PS", "0|1:4000"] and body["rs2"][8:] == ["GT:PS", "1|0:4000"]
    assert body["rs3"][3:5] == ["CAT", "C,CATAT"] and body["rs3"][9] == "2/2"       # the insertion is PharmCAT's second allele
    assert "rs4" not in body                                                          # too few reads: left missing
    assert body["rs5"][9] == "0/0"                                                    # absent from a whole genome: reference
    assert (stats["called"], stats["reference_assumed"], stats["low_quality"], stats["written"]) == (3, 1, 1, 4)


def test_pharmcat_input_does_not_invent_reference_calls_for_partial_data(tmp_path, pc_positions):
    stats, body = pc_rows(tmp_path, [("chr1", 5000, "A", "G", "0/1:15,15:30:99")], whole=False)
    assert set(body) == {"rs1"} and stats["not_covered"] == 4


def test_pharmcat_input_treats_a_reference_mismatch_as_missing(tmp_path, pc_positions):
    stats, body = pc_rows(tmp_path, [("chr1", 5000, "C", "G", "0/1:15,15:30:99")])
    assert "rs1" not in body and stats["mismatched"] == 1


def test_pharmcat_input_passes_through_an_unlisted_snv(tmp_path, pc_positions):
    stats, body = pc_rows(tmp_path, [("chr1", 5000, "A", "T", "0/1:15,15:30:99")])
    assert body["rs1"][4] == "G,T" and body["rs1"][9] == "0/2" and stats["novel_alleles"] == 1


def test_pharmcat_gvcf_uses_blocks_and_skips_uncovered_positions(tmp_path, pc_positions):
    block = ("chr1", 4000, "N", "<NON_REF>", "0/0:30:90:28", ".", "GT:DP:GQ:MIN_DP", "END=6100")
    stats, body = pc_rows(tmp_path, [block])
    assert set(body) == {"rs1", "rs2", "rs3"} and stats["reference_read"] == 3 and stats["not_covered"] == 2


# --------------------------------------------------------------------------------------------- quality, ROH

def test_whole_genome_is_decided_by_snv_count(tmp_path):
    d = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")])
    assert not seqqc.is_whole_genome(d)
    q = seqqc.summary(d)
    assert q["kind"] == "wgs" and not q["whole_genome"] and "exome" in q["flags"][0]
    d.stats["snvs"] = 3_500_000
    assert seqqc.is_whole_genome(d)


def test_roh_from_heterozygous_density(tmp_path):
    from array import array
    d = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")])
    hom = array("I", range(1_000_000, 30_000_000, 2_000))
    het = array("I", [p for p in range(1_000_500, 30_000_000, 1_500) if not 10_000_000 <= p < 16_000_000])
    d.hom_pos, d.het_pos = {"1": hom}, {"1": het}
    roh = seqqc.find_roh(d)
    assert len(roh) == 1 and roh.iloc[0].start == 10_000_001 and roh.iloc[0].end == 16_000_000 and roh.iloc[0].length_mb == 6.0
    d.het_pos = {"1": array("I", range(1_000_500, 30_000_000, 1_500))}
    assert seqqc.find_roh(d).empty


def test_quality_summary_flags_an_implausible_genome(tmp_path):
    d = load(tmp_path, [("chr1", 100, "A", "G", "0/1:15,15:30:99")])
    d.stats.update(snvs=3_500_000, ti_tv=1.4, het_hom_ratio=3.6, mean_depth_at_variants=9.0, passing_alleles=80, failed_alleles=20)
    flags = " ".join(seqqc.summary(d)["flags"])
    assert "Transition/transversion" in flags and "contamination" in flags and "failed the quality gates" in flags and "9x" in flags


# -------------------------------------------------------------------------------------- structural variants

SV_HEADER = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n" + COLS


def sv_file(tmp_path, lines, name="sv.vcf"):
    p = tmp_path / name
    p.write_text(SV_HEADER + "".join(l + "\n" for l in lines))
    return p


def test_structural_variants_are_read_by_type_size_and_zygosity(tmp_path):
    p = sv_file(tmp_path, ["chr17\t43090000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=43100000;SVLEN=-10000\tGT\t0/1",
                           "chr1\t109690000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=109709000\tGT\t1/1",
                           "chr3\t5000000\t.\tN\t<DEL>\t.\tLowQual\tSVTYPE=DEL;END=5000900\tGT\t0/1",
                           "chr9\t2000000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2000500\tGT\t0/0",          # not carried
                           "chr4\t1000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1020\tGT\t0/1",               # under 50 bp
                           "chr5\t1000000\t.\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=300\tGT\t0/1"])
    events, meta = sv.parse(p)
    assert meta["build"] == 38 and meta["records"] == 6
    assert [(e.kind, e.length, e.passed, e.genotype) for e in events] == [
        ("DEL", 10000, True, "heterozygous"), ("DEL", 19000, True, "homozygous"), ("DEL", 900, False, "heterozygous"), ("INS", 300, True, "heterozygous")]


def test_copy_number_records_become_deletions_and_duplications(tmp_path):
    p = sv_file(tmp_path, ["chr22\t42120000\t.\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=42135000\tGT:CN\t./1:3",
                           "chr2\t1000000\t.\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=1100000\tGT:CN\t0/0:2",
                           "chrX\t31000000\t.\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=31200000\tGT:CN\t1:0",
                           "chrX\t40000000\t.\tN\t<CNV>\t.\tPASS\tSVTYPE=CNV;END=40200000\tGT:CN\t1:1"])
    male = [(e.chrom, e.kind, e.genotype, e.copy_number) for e in sv.parse(p, sex="Male")[0]]
    assert male == [("22", "DUP", "", 3), ("X", "DEL", "hemizygous", 0)]      # one copy of X is normal for a male
    female = [(e.chrom, e.kind) for e in sv.parse(p, sex="Female")[0]]
    assert ("X", "DEL") in female and len(female) == 3


def test_structural_findings_follow_dosage_sensitivity(tmp_path):
    p = sv_file(tmp_path, ["chr17\t43090000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=43100000\tGT\t0/1",       # BRCA1: haploinsufficient
                           "chr7\t117500000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=117520000\tGT\t0/1",      # CFTR: recessive
                           "chr1\t109690000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=109709000\tGT\t1/1",      # GSTM1: common, benign
                           "chr17\t14194598\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=15519638\tGT\t0/1",       # the whole CMT1A region
                           "chr17\t14200000\t.\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=14260000\tGT\t0/1",       # a sliver of it
                           "chr3\t1000000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1000400\tGT\t0/1"])         # nowhere of note
    res = sv.analyze([p], 38, "Female")
    f = res["findings"].set_index("gene")
    assert "haploinsufficient" in f.loc["BRCA1", "meaning"] and f.loc["BRCA1", "rank"] == 0
    assert "carrier" in f.loc["CFTR", "meaning"]
    assert "not a disease finding" in f.loc["GSTM1", "meaning"]
    region = [g for g in f.index if g.startswith("17p12")]
    assert len(region) == 1 and f.loc[region[0], "rank"] == 0 and "duplication" in f.loc[region[0], "meaning"]
    assert res["stats"]["passed"] == 6 and res["section"]["id"] == "sv" and len(res["section"]["rows"]) == len(f)
    notable, _ = sv.glance_bullets(res)
    assert any("BRCA1" in n for n in notable)


def test_structural_file_on_another_build_is_skipped(tmp_path):
    p = tmp_path / "old.vcf"
    p.write_text("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n" + COLS + "17\t41200000\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=41250000\tGT\t0/1\n")
    res = sv.analyze([p], 38)
    assert res["findings"].empty and "GRCh37" in res["notes"][0]
    assert sv.glance_bullets(res)[1]


# ------------------------------------------------------------------------------------------- lookup, end to end

def test_fetch_returns_the_records_over_a_position(tmp_path):
    p = vcf(tmp_path, [("chr1", 4000, "N", "<NON_REF>", "0/0:30:90:28", ".", "GT:DP:GQ:MIN_DP", "END=4999"),
                       ("chr1", 5000, "A", "G", "0/1:15,15:30:99"), ("chr1", 9000, "C", "T", "1/1:0,30:30:99")])
    (r,) = fetch(p, "1", 5000, 5000)
    assert (r["ref"], r["alts"], r["genotype"], r["allele_depths"]) == ("A", ["G"], "0/1", "15,15")
    (b,) = fetch(p, "1", 4500, 4500)
    assert b["reference_block"] and b["end"] == 4999
    assert fetch(p, "1", 7000, 7000) == [] and fetch(p, "2", 5000, 5000) == []


def test_cli_runs_a_sequencing_vcf_end_to_end(tmp_path):
    from gnomancer.cli import main
    sites = load_sites(38)
    apoe, lct = sites["rs429358"], sites["rs4988235"]
    p = vcf(tmp_path, [(f"chr{lct.chrom}", lct.pos, lct.ref, "A", "0/1:15,15:30:99"), (f"chr{apoe.chrom}", apoe.pos, apoe.ref, "C", "0/1:15,15:30:99")],
            header=HEADER38, name="genome.vcf")
    out = tmp_path / "out"
    main([str(p), "--out", str(out), "--no-setup", "--no-clinvar", "--no-pharmcat", "--no-y", "--no-mt", "--no-ancestry"])
    md = (out / "REPORT.md").read_text()
    assert md.startswith("# Genome analysis report") and "one e4 allele" in md and "persistent" in md
    assert "not a whole genome" in md                      # two variants: treated as partial data, absent sites unknown
    assert "e3/e4" not in md                               # APOE also needs rs7412, which partial data cannot show as reference
    panel = pd.read_csv(out / "marker-panel.tsv", sep="\t")
    assert (panel.genotype != "no call").sum() == 2
    import json
    s = json.loads((out / "summary.json").read_text())
    assert s["kind"] == "partial" and s["build"] == 38
    assert "<h1><span>gnomancer</span> · genome report</h1>" in (out / "REPORT.html").read_text()


def test_mitochondrial_heteroplasmy_vcf_adds_low_level_variants(tmp_path):
    from gnomancer.seq import merge_mito_vcf
    d = load(tmp_path, [("chrM", 263, "A", "G", "1/1:0,2000:2000:99")])
    mito = tmp_path / "mito.vcf"
    mito.write_text("##fileformat=VCFv4.2\n" + COLS +
                    "chrM\t263\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t1/1:2,1998:2000:0.999\n"          # already in the main VCF
                    "chrM\t3243\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:9600,400:10000:0.04\n"      # low-level, new
                    "chrM\t8993\t.\tT\tG\t.\tstrand_bias\t.\tGT:AD:DP:AF\t0/1:9000,1000:10000:0.1\n" # filtered by the caller
                    "chrM\t11778\t.\tG\tA\t.\tPASS\t.\tGT:AD:DP\t0/1:1500,500:2000\n"                # share from allele depths
                    "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:AD:DP:AF\t0/1:5,5:10:0.5\n")
    assert merge_mito_vcf(d, mito) == 2
    assert d.calls[key("MT", 3243)][0].ab == 0.04 and d.calls[key("MT", 11778)][0].ab == 0.25
    assert key("MT", 8993) not in d.calls and len(d.calls[key("MT", 263)]) == 1


def test_family_mode_declines_sequencing_vcfs_with_a_clear_message(tmp_path):
    from gnomancer.family import run
    a = vcf(tmp_path, [("chr1", 1000, "A", "G", "0/1:10,12:22:99")], name="mum.vcf")
    b = vcf(tmp_path, [("chr1", 1000, "A", "G", "0/1:10,12:22:99")], name="me.vcf")
    with pytest.raises(SystemExit, match="does not read sequencing VCFs"):
        run([str(a), str(b), "--out", str(tmp_path / "fam")])


def test_a_refcall_record_counts_as_an_observed_reference_in_a_gvcf(tmp_path):
    """DeepVariant writes a candidate it rejected as 0/0 with FILTER=RefCall. In a gVCF that is an observed reference
    genotype when depth and quality are adequate, and still a no-call when they are not."""
    block = ("chr1", 4000, "N", "<*>", "0/0:30:50:30", ".", "GT:DP:GQ:MIN_DP", "END=4999")
    weak = Site("rsWeak", "panel", "1", 6000, "C", ("T",))
    d = load(tmp_path, [block, ("chr1", 5000, "A", "G,<*>", "0/0:31,2,0:33:45", "RefCall"), ("chr1", 6000, "C", "T,<*>", "0/0:30,3,0:33:9", "RefCall")], [SNV, weak])
    assert gt(d, SNV) == ("AA", READ)
    assert gt(d, weak) == (None, "covered too thinly to call")
