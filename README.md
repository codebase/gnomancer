# gnomancer

**https://gnomancer.com**

**A private, local deep dive into your raw DNA file.** One command turns an AncestryDNA, 23andMe, MyHeritage,
FamilyTreeDNA, LivingDNA or VCF export into a readable HTML report: file quality, Y and mitochondrial lineages,
traits, carrier sites and pharmacogenomics, each with plain-language context and honest limits.

Nothing is uploaded. The only network traffic is the one-time download of public tools and reference files.

```bash
curl -fsSL https://gnomancer.com/install.sh | bash     # macOS / Linux / Windows via WSL; or clone and run ./setup.sh
gnomancer ~/Downloads/AncestryDNA.zip
open AncestryDNA_analysis/REPORT.html
```

Not sure how to get the file? [gnomancer.com/get-your-data](https://gnomancer.com/get-your-data/) has the steps for every vendor.
Full install notes, including Windows: [gnomancer.com/install](https://gnomancer.com/install/).

Using Claude? The repo doubles as a Claude plugin: `skills/gnomancer/SKILL.md` teaches the agent to install, run and summarize.

![demo report](docs/demo-report.png)

## What you get

`<input>_analysis/`

| File | Contents |
|---|---|
| `REPORT.html` | The report: self-contained, light/dark, printable, searchable marker tables |
| `REPORT.md` | The same content as Markdown |
| `summary.json` | Everything machine-readable |
| `chromosome-qc.tsv`, `roh.tsv`, `marker-panel.tsv` | Raw tables |
| `clinvar-hits.tsv`, `pharmcat-calls.tsv` | ClinVar and PharmCAT results |
| `yhaplo/`, `haplogrep3/`, `pharmcat/` | Raw tool outputs |

Report sections:

0. **In plain language / At a glance** – a narrative summary (yours or Claude's) and an automatic digest: who the data says you are, what is worth knowing, what is reassuring, sensible next steps
1. **Overview** – vendor, build, sex, lineages, blood group, quality flags
2. **Quality** – no-call rate, heterozygosity per chromosome, independent sex inference, contamination flags
3. **Runs of homozygosity** – parental relatedness and founder-population signal
4. **Lineages** – Y haplogroup via [yhaplo](https://github.com/23andMe/yhaplo) with a second opinion from [Yleaf](https://github.com/genid/Yleaf) on current YFull/FTDNA/ISOGG trees; mtDNA via [Haplogrep3](https://github.com/genepi/haplogrep3); each with population context
5. **Multi-site readings** – ABO blood group, lactase persistence, TAS2R38 taster haplotype, APOE (collapsed by default)
6. **Marker panel** – ~110 curated sites: ancestry-informative, traits, health, carrier, pharmacogenomics
7. **ClinVar scan** – every typed site checked against [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/); pathogenic alleles listed with a plausibility score (strand-ambiguous and homozygous-rare calls are flagged as probable array artifacts)
8. **PharmCAT** – star-allele diplotypes and CPIC phenotypes for 18 pharmacogenes via [PharmCAT](https://github.com/PharmGKB/PharmCAT), with per-gene coverage so you can see how much of each gene the chip actually typed
9. **Ancestry estimate** (`--ancestry`) – 1000 Genomes population classification via [ezancestry](https://github.com/arvkevi/ezancestry)
10. **Limits** – what a genotyping array cannot tell you

Health and carrier sections are collapsed by default so that people who would rather not see, say, their APOE
genotype are not shown it by accident.

## Vendors and quirks

Parsing is delegated to [`snps`](https://github.com/apriha/snps), which detects the vendor from the header and
normalizes everything. Quirks that are handled:

| Vendor | Quirk | Handling |
|---|---|---|
| AncestryDNA | CRLF line endings; chromosomes coded 23=X, 24=Y, 25=PAR, 26=MT; two allele columns; no-call `0`; indels `I`/`D`; diploid-looking `GG` on Y and MT | normalized to X/Y/PAR/MT, haploid single letters, NaN no-calls |
| AncestryDNA | V1 (~700k sites, 2012–2016) and V2 (~665k, 2016+) type different sites | report marks `not on array` per site |
| 23andMe | single genotype column; `--` no-call; haploid `A`; internal `i3000001` IDs; v3/v4/v5 chips differ (v5 dropped many PGx sites) | handled; missing sites reported |
| MyHeritage, FamilyTreeDNA | quoted CSV `RSID,CHROMOSOME,POSITION,RESULT` | handled |
| FamilyTreeDNA (old) | build 36 exports | remapped to 37; mtDNA call skipped (build 36 MT is not rCRS) |
| LivingDNA, Genes for Good, tellmeGen, Circle | tab/CSV variants | handled |
| Build 38 exports, VCF | | remapped to GRCh37 (chain file downloaded once) |
| Any | occasional probes reported on the reverse strand | panel auto-flips unambiguous SNPs and flags A/T and C/G ones |
| Any | MT numbering is rCRS (NC_012920) on builds 37/38 | diffed against the bundled rCRS |

## The marker panel

`data/markers.tsv` is a plain TSV anyone can extend:

```
rsid  category  gene  label  ref  effect  interp0  interp1  interp2  note
```

`ref`/`effect` are the two alleles on the **GRCh37 forward strand**. `effect` is the allele that gets counted and
`interp0/1/2` is shown for 0/1/2 copies. Categories: `ancestry`, `trait`, `health`, `carrier`, `pgx`.
See [CONTRIBUTING.md](CONTRIBUTING.md) for the sourcing rules.

Haplogroup context lives in `data/haplogroup_context.tsv`, matched by prefix so a call like `R-M269 / R1b1a2a1a2`
shows the chain R → R1b → R-M269.

## Options

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

## Tools used under the hood

| Tool | What it contributes | Source |
|---|---|---|
| `snps` | vendor detection, normalization, build remapping | github.com/apriha/snps |
| yhaplo | Y haplogroup (ISOGG 2016 tree) | github.com/23andMe/yhaplo |
| Haplogrep3 + PhyloTree 17 | mtDNA haplogroup | github.com/genepi/haplogrep3 |
| ClinVar | clinically classified variants at every typed site | ncbi.nlm.nih.gov/clinvar |
| PharmCAT | CPIC/PharmGKB star-allele calling | github.com/PharmGKB/PharmCAT |
| ezancestry | 1000 Genomes population classification | github.com/arvkevi/ezancestry |
| Yleaf | second-opinion Y haplogroup on current YFull / FTDNA / ISOGG trees | github.com/genid/Yleaf |
| rCRS (NC_012920.1) | mitochondrial reference | NCBI |

Evaluated and not used (yet): PLINK 2 (nothing an array file needs that is not already here), ADMIXTURE (needs a multi-GB reference panel),
pgsc_calc polygenic scores (arrays need imputation first), SNPedia (non-commercial license, no bulk use).

## Try it without your own data

```bash
.venv/bin/python examples/make_demo.py demo      # synthetic AncestryDNA / 23andMe / MyHeritage files
.venv/bin/gnomancer demo/demo_AncestryDNA.txt
```

## Limits, stated plainly

- A genotyping array covers roughly 0.02% of the genome at pre-chosen common sites. It is not a sequence.
- Negative at a site ≠ gene is normal. BRCA1/2, CFTR and similar are covered at a handful of sites only.
- Positive at a rare site is often wrong: 40% of variants from consumer raw data were false positives when a clinical
  lab re-tested them ([Tandy-Connor et al. 2018](https://doi.org/10.1038/gim.2018.38)). The ClinVar section scores each hit for plausibility for this reason.
- CYP2D6 copy number, full HLA typing and rare variants are invisible.
- Y and mtDNA resolution stops at the few thousand / few hundred sites typed.
- Admixture percentages need a reference panel; the vendor's own estimate is the easy source.
- **This is not medical advice.** Confirm anything that matters with a clinical-grade test and a clinician or genetic counselor.

## Website

`site/` is the landing page at https://gnomancer.com, served as a Cloudflare Worker with static assets (`wrangler.toml`).
Deploy with `npx wrangler deploy` (Node 22+, `wrangler login` once). It is static: no cookies, no analytics.

## Development

```bash
.venv/bin/pip install -e ".[dev]"
.venv/bin/pytest
```

CI runs the tests and builds a demo report on every push. Real raw-data files and `*_analysis/` folders are
git-ignored; never commit them.

## Related tools

Local: [`apriha/lineage`](https://github.com/apriha/lineage) (shared segments between kits),
[PLINK 2](https://github.com/chrchang/plink-ng), ADMIXTURE + 1000 Genomes, [SNPedia](https://www.snpedia.com).
Web (uploads your data): GEDmatch, James Lick's mthap, Morley's Y-SNP predictor, Promethease, Genetic Genie.

## License

MIT. Haplogroup tools and references have their own licenses: yhaplo (23andMe), Haplogrep3 (GENEPI Innsbruck),
PhyloTree, ISOGG.
