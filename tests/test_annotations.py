import gzip
from pathlib import Path

import pandas as pd

from gnomancer import clinvar, pharmcat

ROOT = Path(__file__).resolve().parent.parent


def _mini_clinvar(tmp_path):
    p = tmp_path / "clinvar.vcf.gz"
    lines = ["##fileformat=VCFv4.1", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             "1\t100\t1\tC\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=GENEA:1;CLNDN=Disease_A;RS=1001;AF_EXAC=0.0002",
             "1\t200\t2\tA\tT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;GENEINFO=GENEB:2;CLNDN=Disease_B;RS=1002",
             "1\t300\t3\tG\tA\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENEC:3;RS=1003",
             "1\t400\t4\tG\tA\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;GENEINFO=GENED:4;CLNDN=Disease_D;RS=1004;AF_TGP=0.12",
             "1\t500\t5\tGA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=practice_guideline;GENEINFO=GENEE:5;RS=1005"]
    with gzip.open(p, "wt") as f:
        f.write("\n".join(lines) + "\n")
    return p


def test_clinvar_scan_flags_and_filters(tmp_path):
    vcf = _mini_clinvar(tmp_path)
    df = pd.DataFrame({"chrom": ["1"] * 5, "pos": [100, 200, 300, 400, 500],
                       "genotype": ["CT", "TT", "AA", "AA", "DI"]}, index=["rs1001", "rs1002", "rs1003", "rs1004", "rs1005"])
    r = clinvar.scan(df, vcf)
    h = r["hits"].set_index("rsid")
    assert r["sites_in_clinvar"] == 5 and r["pathogenic_sites_typed"] == 4
    assert set(h.index) == {"rs1001", "rs1002", "rs1004"}          # benign and indel rows excluded
    assert h.loc["rs1001"].confidence == "high" and h.loc["rs1001"].stars == 3
    assert h.loc["rs1002"].confidence == "low" and "strand-ambiguous" in h.loc["rs1002"]["flags"]
    assert h.loc["rs1004"].confidence == "medium" and "common" in h.loc["rs1004"]["flags"]
    assert "3 site(s)" in clinvar.interpret(r)


def test_clinvar_missing_file_is_graceful(tmp_path):
    r = clinvar.scan(pd.DataFrame(columns=["chrom", "pos", "genotype"]), tmp_path / "nope.vcf.gz")
    assert not r["available"] and "gnomancer setup" in clinvar.interpret(r)


def test_pharmcat_vcf_builder(tmp_path, monkeypatch):
    pos = tmp_path / "positions.vcf"
    pos.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n"
                   "chr10\t94942290\trs1799853\tC\tT\t.\tPASS\tPX=CYP2C9\tGT\t0/0\n"
                   "chr10\t94981296\trs1057910\tA\tC\t.\tPASS\tPX=CYP2C9\tGT\t0/0\n"
                   "chr10\t1\trs999\tA\tG\t.\tPASS\tPX=GENEX\tGT\t0/0\n"
                   "chr10\t2\trs998\tAT\tA\t.\tPASS\tPX=GENEX\tGT\t0/0\n")
    monkeypatch.setattr(pharmcat, "POSITIONS", pos)
    df = pd.DataFrame({"chrom": ["10"] * 4, "pos": [1, 2, 3, 4], "genotype": ["CT", "AT", "TC", "DI"]},
                      index=["rs1799853", "rs1057910", "rs999", "rs998"])
    stats = pharmcat.build_vcf(df, tmp_path / "s.vcf")
    body = [l for l in (tmp_path / "s.vcf").read_text().splitlines() if not l.startswith("#")]
    assert stats["written"] == 2 and stats["mismatched"] == 1 and stats["indels_skipped"] == 1
    assert body[0].endswith("0/1")                      # CT -> het
    assert body[1].split("\t")[2] == "rs999" and body[1].endswith("0/1")   # TC flipped to AG -> het


def test_yleaf_plink_and_parse(tmp_path):
    from gnomancer import yleaf
    df = pd.DataFrame({"chrom": ["Y", "Y", "Y", "1"], "pos": [100, 200, 300, 5], "genotype": ["G", "AG", "TT", "AA"]},
                      index=["rs1", "rs2", "rs3", "rs4"])
    n = yleaf.write_plink(df, tmp_path / "s")
    assert n == 2 and (tmp_path / "s.map").read_text().count("\n") == 2
    assert (tmp_path / "s.ped").read_text().strip().endswith("G\tG\tT\tT")
    hg = tmp_path / "c.hg"
    hg.write_text("Sample_name\tyfull_Hg\tyfull_Hg_marker\tyfull_Valid_markers\tyfull_QC\tisogg_Hg\tisogg_Hg_marker\tisogg_Valid_markers\tisogg_QC\n"
                  "s\tR-L21*\tL21\t900\t1.0\tR1b1a2a1a2c\tL21\t850\t0.9\n")
    t = yleaf.parse_combined(hg)
    assert t["yfull"]["haplogroup"] == "R-L21*" and t["isogg"]["qc"] == "0.9" and "ftdna" not in t
