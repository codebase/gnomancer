# Security

gnomancer runs entirely on the user's machine and never transmits the genotype file. The attack surface is the
local pipeline (parsing an untrusted file) and the optional MCP server (stdio, local only, no network listener).

**Reporting.** Email security reports to the maintainers via a private GitHub security advisory on this repository
(Security tab → Report a vulnerability). Please do not open public issues for vulnerabilities. We aim to acknowledge
within 72 hours.

**Scope of interest:** path traversal or code execution via crafted raw-data files, the installer script, the MCP
tools (e.g. `analyze` or `lookup` reading files outside the user's intent), and the downloaded third-party tools'
verification (we pin versions but do not yet verify checksums; see roadmap).

**Out of scope:** the biology (interpretation accuracy is a correctness issue, filed as a normal issue), and the
third-party tools themselves (report those upstream: yhaplo, Haplogrep3, PharmCAT, Yleaf, snps).
