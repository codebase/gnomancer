# Referrals and roadmap notes

Working notes, not commitments. Last checked 2026-11-02.

## Referral / affiliate programs that fit the report

The report points people to a next step in four places (Y sequencing, confirmatory testing, a genetic counselor, Lp(a) blood test). These are the programs that could sit behind a future `/next-steps` page; none is linked today.

| Report says | Vendor | Program | Terms found | Fit |
|---|---|---|---|---|
| "Y sequencing resolves deeper than an array" | FamilyTreeDNA (Big Y-700, mtFull Sequence) | CJ Affiliate ([program page](https://cms.familytreedna.com/?p=29213), [FlexOffers listing](https://www.flexoffers.com/affiliate-programs/familytreedna-affiliate-program/)) | ~4% per sale, 7-day cookie, $50 payout threshold; Big Y ≈ $450 so ≈ $18/sale | Best fit: genealogy, no medical claim |
| "Arrays miss rare variants" | Sequencing.com (30x WGS) | [FlexOffers](https://www.flexoffers.com/affiliate-programs/sequencing-com-affiliate-program/) | 4% per sale, 45-day cookie | Good: honest fix for array limits |
| same | Nebula Genomics / DNA Complete | [Awin](https://ui.awin.com/merchant-profile/65972) | 1-day cookie, rate not shown | OK |
| "Confirm with a clinical test / talk to a counselor" | Genome Medical | [provider referral](https://www.genomemedical.com/provider-referral/) | partnership channel, no consumer payout found | Link anyway: it is the trust feature |
| "Try a second vendor" (visitors without data) | MyHeritage, CircleDNA | own program / [Commission Factory](https://www.commissionfactory.com/advertiser-directory/circledna-affiliate-program/85028) | exist, not verified | Low priority |
| — | 23andMe | [FlexOffers](https://www.flexoffers.com/affiliate-programs/23andme-com-affiliate-program) | "not currently offering" | Skip |
| — | Dante Labs | — | reportedly out of business since 2024 | Skip |
| — | consumer PGx tests (ClarityX etc.) | — | nothing solid found; most PGx is physician-ordered | Skip |
| "Lp(a) / ferritin / HbA1c is a cheap blood test" | DTC lab marketplaces | — | not researched yet | Maybe |

### Rules for using any of them
1. Genealogy upsells only next to genealogy content. Never an affiliate link beside a disease-risk row.
2. If this is ever done, the report links to **one page**, `gnomancer.com/next-steps` (not built yet); affiliate URLs live there, not in the tool. Partners can change without a release, clicks are measurable, the open-source report stays clean.
3. Plain-words disclosure on that page. Keep the genetic-counselor option on it even though it earns nothing.
4. Expect pocket money: at 4%, a thousand runs a month is a handful of sales. It is the only monetization that does not contradict "nothing leaves your machine".

### Sign-up (needs tax and payout details, so the maintainer does this)
- FamilyTreeDNA: apply via CJ Affiliate, program listed on their page above.
- Sequencing.com: FlexOffers.
- Nebula / DNA Complete: Awin.

## Distribution and discoverability (for people without a GitHub login)

Nobody needs a GitHub account to use gnomancer: the repo is public, releases are plain downloads, and the installer is
`curl -fsSL https://gnomancer.com/install.sh | bash`. The gaps are packaging and being found.

**Packaging:**
- ☑ `install.sh` on gnomancer.com and the `/install/` page (macOS, Linux, Windows via WSL).
- ☑ PyPI (`pipx install gnomancer`), trusted publishing from GitHub releases.
- ☑ Homebrew tap (`brew install codebase/tap/gnomancer`, python@3.13 + openjdk); unverifiable on macOS 26.2 until Apple's next point release, see development.md.
- ☑ Claude skill and plugin (`skills/gnomancer/SKILL.md`, `.claude-plugin/plugin.json`).
- ☑ MCP server on the official registry plus Glama, Smithery, mcp.so, PulseMCP (docs/mcp-registry.md).
- ⏸ Native Windows and a desktop app: parked (see roadmap item 7).

**Where people look (announce in this order):**
- r/AncestryDNA, r/23andme, r/Genealogy, r/promethease, r/SNPedia: the "what can I do with my raw data" threads are weekly.
- ISOGG wiki "Autosomal DNA tools" and the SNPedia "tools" page: standing lists people actually consult.
- Genetic-genealogy blogs that review tools: DNAeXplained (Roberta Estes), The DNA Geek, Your DNA Guide, Genealogy Explained.
- Hacker News "Show HN" and Product Hunt.
- A 10-minute YouTube walkthrough: download the file, run, read the report. Search traffic for "what to do with raw DNA data" is large.

**SEO pages on the site** (each is a question people type): `/get-your-data/` ☑, `/install/` ☑, one page per haplogroup ☑ (67 pages), "what is a raw DNA file" ☐, "is my 23andMe data safe to upload to X" ☐.

## Community

- ☑ Discord is the front door for people without a GitHub login: https://discord.gg/aGjqzqErQH (verification, rules screening, channels per topic and haplogroup). Linked from the site footer and README.
- GitHub Issues stay for bugs and developer talk.
- House rule, posted in the server: share haplogroups, traits and questions; never paste raw files, ClinVar rows or health tables. Genetic data posted publicly is public forever.
- A hosted Discourse (email sign-up, Google-searchable threads) is the step after Discord if the community outgrows it; about $100/month hosted.

## Roadmap (prioritized 2026-11-02)

Status key: ☐ not started · ◐ in progress · ☑ done

### Make the report better (private by construction)
1. ☑ **Population frequencies on every panel row.** "About 31% of Europeans share this genotype." 1000 Genomes super-population frequencies baked into `marker_frequencies.tsv`; genotype frequency under Hardy-Weinberg for the user's ancestry group (or EUR by default). Half a day.
2. ☑ **ClinVar reclassification watch.** `gnomancer watch <report>` re-downloads ClinVar when NCBI publishes a new release, rescans the user's typed sites, and reports only what changed (new pathogenic hits, upgrades, downgrades). Snapshot stored with the report. Half a day.
3. ☑ **Pharmacist wallet card.** One-page PDF of the PharmCAT table with the "confirm clinically" footer. An hour.
4. ☑ **Family mode** (`gnomancer family`, MCP tool `family`): shared segments via lineage, relationship estimate, shared lineages, trio inheritance and Mendelian check. Synthetic family generator in `examples/make_family.py`.
5. ☐ **Polygenic risk scores** (PGS Catalog) for a few well-validated traits as percentiles vs 1000 Genomes, with the imputation caveat. A week; more medical-adjacent.

### Reach more people

The assistant path (MCP server + skill) is the front door; `analyze` now auto-runs setup on first use so there is no manual step at all.
6. ⏸ **Run in the browser** (Pyodide). Parked 2026-11-02: a web page asking for a genome contradicts the product's one promise, however it works under the hood. Revisit only as a no-file-picker demo on the synthetic sample.
7. ⏸ **Native desktop app**. Parked 2026-11-02: needs paid signing certificates (Apple $99/yr, Windows $200+/yr) to avoid scary warnings, and the audience it serves is reached better through the MCP server and skill, which auto-install. Revisit if non-assistant users ask for it.
8. ☑ **Homebrew tap** (`brew install codebase/tap/gnomancer`; formula in github.com/codebase/homebrew-tap, bump url/sha256 per release) that pulls python and openjdk.
10. ☑ **Haplogroup story pages**: 67 pages under gnomancer.com/haplogroups/ generated from `haplogroup_context.tsv` by `site/build_haplogroups.py` (re-run after editing the table); report lineage cards link to them. No map yet.

### Polish and trust
11. ☑ **Checksum verification** for downloaded tools (SHA-256 pins; ClinVar against NCBI's MD5 sidecar).
12. ☑ **Light-mode and print-to-PDF pass** on the report (checked 2026-11-02; print CSS collapses controls and expands the health sections).
13. ☐ **Walkthrough video** for the site and YouTube. Needs the maintainer on mic, or a script.

### Done
- ☑ Core pipeline, HTML/Markdown reports, at-a-glance digest, narrative slot (`--summary`, `--ai`)
- ☑ ClinVar scan with plausibility scoring; PharmCAT; yhaplo + Yleaf; Haplogrep3; ezancestry
- ☑ gnomancer.com (home, get-your-data, install, installer script), PyPI, MCP server + official registry, Glama, mcp.so and awesome-list submissions, Claude skill/plugin
- ☑ Logo, brand accent, README rewrite

### Waiting on the maintainer
- ☑ Discord server: https://discord.gg/aGjqzqErQH (Medium verification, Community mode, rules screening, 12 channels)
- Affiliate sign-ups (FamilyTreeDNA via CJ, Sequencing.com via FlexOffers) if `/next-steps` should earn
- ISOGG / SNPedia / DNA Painter accounts and the Reddit posts
