# Command-line options

```
gnomancer FILE [--out DIR] [--markers TSV] [--no-y] [--no-yleaf] [--no-mt] [--no-clinvar] [--no-pharmcat] [--ancestry]
          [--summary summary.md] [--ai]
```

**Plain-language summary.** The report opens with an automatic "at a glance" block built from the results. For a narrative on
top of that, every run also writes `summary-prompt.md`: paste it to whatever AI assistant you use (Claude, ChatGPT, a local
model), save the answer as `summary.md`, and re-run with `--summary summary.md`. Or pass `--ai` to have Claude write it
directly (`pip install anthropic`, credentials via `ANTHROPIC_API_KEY` or `ant auth login`). `--ai` is the one path that sends
anything off your machine: the structured results in `summary-prompt.md`, never the raw genotype file.

`--ancestry` runs [ezancestry](https://github.com/arvkevi/ezancestry) for a 1000 Genomes population classification.
ezancestry and Yleaf need Python ≤3.11; `gnomancer setup` builds a side venv (`~/.gnomancer/venv-ancestry`) for both automatically when
`python3.11` is on the PATH, and the CLI uses it when present. First runs download reference data (ezancestry models;
Yleaf's hg19 reference, ~900 MB).

## Try it without your own data

```bash
python examples/make_demo.py demo      # synthetic AncestryDNA / 23andMe / MyHeritage files (from a repo checkout)
gnomancer demo/demo_AncestryDNA.txt
```

## Pharmacist card

Every run also writes `pgx-card.pdf` (and `pgx-card.html`): a one-page summary of the PharmCAT diplotypes, phenotypes and the drugs each gene matters for, with the clotting and abacavir single-site results and a footer that says it is array-based and provisional. Print it or show it on a phone; the report links to it from the PharmCAT section.

## `gnomancer watch REPORT_DIR`: ClinVar reclassification check

ClinVar reclassifies variants every month. `watch` keeps a snapshot of every ClinVar-classified allele you carry
(any significance) next to your report (`clinvar-snapshot.json`, written on every normal run), checks whether NCBI
has published a newer ClinVar release, downloads it if so, rescans the same raw file, and prints only what changed:

- newly pathogenic / likely pathogenic at a site you carry (new record or reclassification)
- downgraded (no longer pathogenic)
- review-status changes on pathogenic sites

The result is also saved as `clinvar-watch.md` in the report folder. Flags: `--force` (rescan even if ClinVar is
unchanged), `--no-download` (use the ClinVar file already on disk). The raw file must still be at the path recorded
in `summary.json`.

Run it monthly. macOS/Linux cron, 9:00 on the 2nd of each month (ClinVar releases on the 1st):

```
0 9 2 * * ~/.gnomancer/venv/bin/gnomancer watch ~/Downloads/AncestryDNA_analysis >> ~/.gnomancer/watch.log 2>&1
```

From an AI assistant with the MCP server connected: ask it to run `clinvar_watch`.
