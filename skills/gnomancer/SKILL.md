---
name: gnomancer
description: Analyze a consumer raw DNA file (AncestryDNA, 23andMe, MyHeritage, FamilyTreeDNA, LivingDNA, VCF) locally with gnomancer and explain the results. Use when the user mentions their DNA raw data, a genotype file, haplogroups, or asks what their 23andMe/Ancestry data says. Never uploads the file anywhere.
---

# gnomancer: local raw-DNA analysis

gnomancer turns a raw genotype export into a report: quality, Y and mtDNA haplogroups, a curated trait/health/
pharmacogenomics panel, a ClinVar scan, PharmCAT star alleles, and a plain-language summary. Everything runs on the
user's machine. Your job is to run it, write the narrative summary it asks for, and answer questions from the output.

## Privacy rules (non-negotiable)
- Never send the raw file, or any genotype rows from it, to any service or paste them into chat. Work from the
  report files gnomancer writes (`summary.json`, `REPORT.md`, `summary-prompt.md`).
- Do not read `clinvar-hits.tsv` or the health sections aloud unless the user asks; the report collapses them on purpose.
- Say plainly that an array is probabilistic and rare-site calls need clinical confirmation. Never diagnose.

## If the gnomancer MCP server is connected
Prefer its tools: `analyze(path)`, then `narrative_prompt()` → write the narrative → `set_narrative(text)`, then `summary()`,
`section(name)`, `lookup(rsid)`. For two or more relatives' files use `family(paths)`. Use `clinvar_findings()` / `pharmacogenomics()` only when asked about health or drugs; `clinvar_watch()` re-checks ClinVar for reclassifications when the user asks whether anything changed.
Otherwise follow the steps below with the CLI.

## Steps
1. **Install if missing.** `command -v gnomancer || ls ~/.gnomancer/venv/bin/gnomancer`. If absent: `pipx install gnomancer && gnomancer setup`
   (or `curl -fsSL https://gnomancer.com/install.sh | bash`; Windows needs WSL, see https://gnomancer.com/install/).
   Prerequisites: Python 3.10+, Java 17+ (Java only for the mtDNA haplogroup and PharmCAT).
2. **Find the file.** Ask for the path if not given. Zips are fine. If the user has no file yet, point them to
   https://gnomancer.com/get-your-data/ for per-vendor download steps.
3. **Run it.** `gnomancer <file> --out <dir>` (add `--ancestry` if the user wants a population estimate; it is slow
   the first time). Expect one to three minutes (the ClinVar scan dominates).
4. **Write the narrative.** Open `<dir>/summary-prompt.md`, follow its instructions exactly (250–400 words, second
   person, no headings, no invented facts, no rsIDs), save your text as `<dir>/summary.md`, then re-run with
   `--summary <dir>/summary.md` so it appears at the top of the report.
5. **Show the result.** Open `<dir>/REPORT.html` for the user and give a five-line recap: lineages, blood group,
   anything in "Worth knowing", anything flagged for clinical confirmation, and the reminder that it is not medical advice.
6. **Answer questions** from `summary.json` (machine-readable: qc, y, mt, combos, clinvar, pharmcat) and `REPORT.md`.
   For "what does site rsXXXX say", grep `marker-panel.tsv`. If a site is not on the chip, say so; do not guess.

## Interpreting
- Y and mtDNA haplogroups: use the population context already in the report; do not over-claim ethnicity from one lineage.
- "effect_count 2/2" on a health row means two copies of the effect allele; the row's reading gives the odds as published.
- ClinVar hits marked `low` confidence are probable array artifacts (homozygous rare or strand-ambiguous); say so.
- PharmCAT "partial" coverage means untyped positions were assumed reference; results are provisional.
