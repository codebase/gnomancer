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

## Community

- **Start with GitHub Discussions** (free, no moderation infra, no accounts with us). Turn on in repo settings.
- **If a forum on gnomancer.com is wanted**, embed Discussions with [giscus](https://giscus.app): it renders as a forum on our domain, uses GitHub login, zero servers. Structure categories by haplogroup (like FTDNA projects) plus "Reading your report" and "Tools".
- House rule to post prominently: share haplogroups, traits and questions; never paste raw files, ClinVar rows or health tables. Genetic data posted publicly is public forever.
- A self-hosted Discourse is the step after that, only if giscus is outgrown.

## Value-add ideas, roughly in order of payoff ÷ effort

1. **Population frequencies on every panel row.** Embed 1000 Genomes / gnomAD allele frequencies per super-population in `markers.tsv`, so each row says "about 31% of Europeans share this genotype". No server, no user data, makes every row meaningful. One afternoon.
2. **ClinVar reclassification watch.** `gnomancer watch` re-downloads ClinVar monthly and diffs classifications at the user's typed sites; emails/prints only what changed. Clinicians do this by hand; nobody offers it to consumers. Half a day.
3. **Run in the browser (Pyodide/WebAssembly).** "Analyze without installing" on gnomancer.com, file never leaves the tab. The biggest adoption unlock and the only honest way to offer a paid tier later. One to two weeks; Java tools (Haplogrep3, PharmCAT) need WASM ports or a reduced browser mode.
4. **MCP server for the report.** `gnomancer serve` exposes the results to Claude Desktop / any MCP client so people can ask "what does my report say about caffeine?" locally. Fits the AI-agent angle, small effort.
5. **Family mode.** Two or more kits: shared segments via `apriha/lineage`, parent/child/sibling confirmation, which parent each lineage came from, inherited-together traits. Medium effort, high emotional value.
6. **Printable pharmacogenomics wallet card.** One-page PDF with the PharmCAT table and "confirm clinically" footer, for the pharmacist. Small effort, very practical.
7. **Privacy scrubber + vendor converter.** Strip health and carrier sites from a raw file before uploading to GEDmatch; convert between vendor formats (snps can). On-brand, unique.
8. **Haplogroup story pages on the site.** One page per major Y/mt haplogroup from `haplogroup_context.tsv`, with a map. SEO, content, and the natural home for per-haplogroup discussion.
9. **Shareable cards for non-sensitive results** (haplogroups, taste, earwax, caffeine) as PNGs for social. Growth loop with health kept out by design.
10. **Polygenic risk scores** from PGS Catalog for a few well-validated traits, reported as percentiles against 1000 Genomes with the imputation caveat. Big "wow", larger effort, more medical-adjacent.
11. **One-click desktop app** (Tauri wrapper around the CLI). Adoption for non-terminal users.
