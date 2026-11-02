# Contributing

Thanks for helping. The two most valuable contributions are **marker rows** and **vendor fixtures**.

## Adding a marker
Edit `gnomancer/data/markers.tsv`. Every row needs:
- `ref` and `effect` on the **GRCh37 forward strand**. Check dbSNP (the "+" strand alleles for GRCh37) or SNPedia's
  orientation field. If the gene is on the minus strand, the literature allele is the complement of what goes here.
- A source for the interpretation text: a replicated GWAS, a CPIC/PharmGKB guideline, ClinVar, or a review. Put it in `note` if it is not obvious.
- Hedged wording. Effect sizes as "~2x odds", never "you will".
- No rows whose only support is a single small study or a commercial trait report.

Run `pytest` afterwards: the test suite checks the table parses, alleles are valid, and there are no duplicates.

## Adding a vendor
Parsing lives in the [`snps`](https://github.com/apriha/snps) library. If a file does not load, first check whether `snps`
supports it; if not, open an issue there and here with the first 20 lines of the file (header only, no genotypes).
Add an anonymised fixture under `tests/fixtures/` generated with `examples/make_demo.py` style synthetic data, never real data.

## Rules
- Never commit a real raw-data file or a real report. `.gitignore` blocks the common names; stay alert anyway.
- Keep the tool offline: the only network calls allowed are downloading public tools and reference files.
- Health wording goes through the "would a genetic counselor wince?" test.
