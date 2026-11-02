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
ezancestry and Yleaf need Python ≤3.11; `setup.sh` builds a side venv (`.venv-ancestry`) for both automatically when
`python3.11` is on the PATH, and the CLI uses it when present. First runs download reference data (ezancestry models;
Yleaf's hg19 reference, ~900 MB).

## Try it without your own data

```bash
.venv/bin/python examples/make_demo.py demo      # synthetic AncestryDNA / 23andMe / MyHeritage files
.venv/bin/gnomancer demo/demo_AncestryDNA.txt
```
