# MCP registries

The server is `gnomancer serve` (stdio). `server.json` at the repo root is the official MCP Registry manifest; the `mcp` package is a core dependency, so `uvx gnomancer serve` works with nothing else installed.

## Official MCP Registry (registry.modelcontextprotocol.io) — published

Ownership is proven two ways: the PyPI README carries the line `mcp-name: com.gnomancer/gnomancer`, and the domain serves an Ed25519 public key at https://gnomancer.com/.well-known/mcp-registry-auth. The private key lives outside the repo.

To publish a new version:
1. Release to PyPI first (see the release checklist in development.md).
2. Bump `version` in `server.json` (both the top-level field and the package entry).
3. `mcp-publisher login http --domain gnomancer.com --private-key <hex seed>` then `mcp-publisher publish` from the repo root.

## Other directories — listed

Glama, Smithery, mcp.so, PulseMCP and the awesome-mcp-servers list were submitted on 2026-11-02. They index the GitHub repo and re-read `server.json`, so no action is needed per release beyond keeping the manifest current.

## Claude Desktop extension (.mcpb) — parked

A one-click bundle needs Python and the tools packaged with it, which is the same work as the parked desktop app. Revisit together.

## Privacy line for every listing
"Runs locally over stdio. Reads only the user's raw file on their machine. Returns derived results and single-site
lookups to the assistant; never the raw file. Health tables are separate tools fetched only on request."
