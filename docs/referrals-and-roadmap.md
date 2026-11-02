# Referrals and roadmap notes

Working notes, not commitments. Last checked 2026-11-02.

## Referral / affiliate programs that fit the report

The report already sends people to a next step in four places. These are the programs behind them.

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
2. The report links to **one page**, `gnomancer.com/next-steps`; affiliate URLs live there, not in the tool. Partners can change without a release, clicks are measurable, the open-source report stays clean.
3. Plain-words disclosure on that page. Keep the genetic-counselor option on it even though it earns nothing.
4. Expect pocket money: at 4%, a thousand runs a month is a handful of sales. It is the only monetization that does not contradict "nothing leaves your machine".

### Sign-up (needs tax and payout details, so the maintainer does this)
- FamilyTreeDNA: apply via CJ Affiliate, program listed on their page above.
- Sequencing.com: FlexOffers.
- Nebula / DNA Complete: Awin.

## Distribution and discoverability (for people without a GitHub login)

Nobody needs a GitHub account to use gnomancer: the repo is public, releases are plain downloads, and the installer is
`curl -fsSL https://gnomancer.com/install.sh | bash`. The gaps are packaging and being found.

**Packaging, in order:**
1. `install.sh` on gnomancer.com (done) and the `/install/` page with macOS, Linux and Windows-via-WSL steps (done).
2. **PyPI**: `pipx install gnomancer`. pyproject is ready; needs a PyPI account and a trusted-publisher GitHub Action. the maintainer creates the account.
3. **Homebrew tap** (`brew install codebase/tap/gnomancer`): a formula that depends on python@3.12 and openjdk; solves the Java problem on macOS.
4. **Native Windows**: PyInstaller build of the CLI plus bundled JRE, or a Tauri desktop app. The desktop app is the real answer for non-terminal users.
5. **Claude skill / plugin** (done: `skills/gnomancer/SKILL.md`, `.claude-plugin/plugin.json`): any Claude Code or Claude Desktop user can say "analyze my DNA file" and the agent installs, runs, and summarizes. Submit to the Claude Code plugin marketplace and the Anthropic skills directory.
6. **MCP server** (done: `gnomancer serve`, `server.json`, `docs/mcp-registry.md`). Registry publication is blocked on PyPI (step 2).

**Where people look (announce in this order):**
- r/AncestryDNA, r/23andme, r/Genealogy, r/promethease, r/SNPedia: the "what can I do with my raw data" threads are weekly.
- ISOGG wiki "Autosomal DNA tools" and the SNPedia "tools" page: standing lists people actually consult.
- Genetic-genealogy blogs that review tools: DNAeXplained (Roberta Estes), The DNA Geek, Your DNA Guide, Genealogy Explained.
- Hacker News "Show HN" and Product Hunt, once the browser or desktop version exists.
- A 10-minute YouTube walkthrough: download the file, run, read the report. Search traffic for "what to do with raw DNA data" is large.

**SEO pages on the site** (each is a question people type): `/get-your-data/` (done), `/install/` (done), one page per major haplogroup, "what is a raw DNA file", "is my 23andMe data safe to upload to X".

## Community

- Most of the audience has no GitHub login, so GitHub Discussions / giscus is the wrong front door for them (fine for developers).
- **Discord** is the pragmatic consumer forum: nearly everyone has an account, free, good moderation tools, channels per haplogroup. Create the server, add an invite link to the site footer.
- A hosted **Discourse** (email sign-up, searchable by Google, threads that last) is the better long-term forum if the community grows; ~$100/month hosted, or self-host.
- Keep GitHub Issues/Discussions for bugs and developer talk.
- House rule to post prominently: share haplogroups, traits and questions; never paste raw files, ClinVar rows or health tables. Genetic data posted publicly is public forever.
- A self-hosted Discourse is the step after that, only if giscus is outgrown.

## Roadmap (prioritized 2026-11-02)

Status key: ☐ not started · ◐ in progress · ☑ done

### Make the report better (private by construction)
1. ◐ **Population frequencies on every panel row.** "About 31% of Europeans share this genotype." 1000 Genomes super-population frequencies baked into `marker_frequencies.tsv`; genotype frequency under Hardy-Weinberg for the user's ancestry group (or EUR by default). Half a day.
2. ☑ **ClinVar reclassification watch.** `gnomancer watch <report>` re-downloads ClinVar when NCBI publishes a new release, rescans the user's typed sites, and reports only what changed (new pathogenic hits, upgrades, downgrades). Snapshot stored with the report. Half a day.
3. ☐ **Pharmacist wallet card.** One-page PDF of the PharmCAT table with the "confirm clinically" footer. An hour.
4. ☐ **Family mode.** Two or more kits: relatedness confirmation, which parent gave which lineage, shared segments via `apriha/lineage`. A few days.
5. ☐ **Polygenic risk scores** (PGS Catalog) for a few well-validated traits as percentiles vs 1000 Genomes, with the imputation caveat. A week; more medical-adjacent.

### Reach more people
6. ☐ **Run in the browser** (Pyodide/WebAssembly). "Analyze without installing" on gnomancer.com, file never leaves the tab. Biggest adoption unlock and the only honest base for a paid tier. One to two weeks; Java tools need a reduced browser mode.
7. ☐ **Native Windows / desktop app** (Tauri wrapper, bundled JRE). A week.
8. ☐ **Homebrew tap** (`brew install codebase/tap/gnomancer`) that pulls python and openjdk. An hour.
10. ☐ **Haplogroup story pages** on the site, one per major lineage, with a map. SEO plus the natural home for community. A day.

### Polish and trust
11. ☐ **Checksum verification** for downloaded tools (SECURITY.md currently admits we don't). An hour.
12. ☐ **Light-mode and print-to-PDF pass** on the report. An hour.
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
