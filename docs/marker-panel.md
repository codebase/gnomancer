# The marker panel

`gnomancer/data/markers.tsv` is a plain TSV anyone can extend:

```
rsid  category  gene  label  ref  effect  interp0  interp1  interp2  note
```

`ref`/`effect` are the two alleles on the **GRCh37 forward strand**. `effect` is the allele that gets counted and
`interp0/1/2` is shown for 0/1/2 copies. Categories: `ancestry`, `trait`, `health`, `carrier`, `pgx`.
See [CONTRIBUTING.md](CONTRIBUTING.md) for the sourcing rules.

Haplogroup context lives in `gnomancer/data/haplogroup_context.tsv`, matched by prefix so a call like `R-M269 / R1b1a2a1a2`
shows the chain R → R1b → R-M269.
