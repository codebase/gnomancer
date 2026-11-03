# The marker panel

`gnomancer/data/markers.tsv` is a plain TSV anyone can extend:

```
rsid  category  gene  label  ref  effect  interp0  interp1  interp2  note
```

`ref`/`effect` are the two alleles on the **GRCh37 forward strand**. `effect` is the allele that gets counted and
`interp0/1/2` is shown for 0/1/2 copies. Categories: `ancestry`, `trait`, `health`, `carrier`, `pgx`.
See [CONTRIBUTING.md](../CONTRIBUTING.md) for the sourcing rules.

Haplogroup context lives in `gnomancer/data/haplogroup_context.tsv`, matched by prefix so a call like `R-M269 / R1b1a2a1a2`
shows the chain R → R1b → R-M269.

## Population frequencies

`gnomancer/data/marker_frequencies.tsv` holds the effect-allele frequency of every panel site in the five 1000 Genomes
phase 3 super-populations (AFR, AMR, EAS, EUR, SAS), fetched from Ensembl's GRCh37 REST API by
`examples/fetch_frequencies.py`. The report turns these into "how common is my genotype" under Hardy-Weinberg
(p², 2pq, q²), shown for the user's super-population when `--ancestry` ran, else EUR, with all five on hover.
Sites without 1000 Genomes data (mostly rare pharmacogenomic alleles) show a dash. Regenerate the table after
adding markers.
