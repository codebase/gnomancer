# Publishing the MCP server to registries

The server is `gnomancer serve` (stdio). `server.json` at the repo root is the official MCP Registry manifest.

## Official MCP Registry (registry.modelcontextprotocol.io)
Prerequisite: the package must be on **PyPI** as `gnomancer`, and its PyPI README must contain the line
`mcp-name: com.gnomancer/gnomancer` so the registry can verify ownership.
1. Publish to PyPI (trusted publisher via GitHub Actions is the clean way; needs a PyPI account owned by the maintainer).
2. Install the publisher CLI: `brew install mcp-publisher` (or download from github.com/modelcontextprotocol/registry/releases).
3. `mcp-publisher login http --domain gnomancer.com --private-key <hex seed>` (ownership is proven by the Ed25519 public key served at https://gnomancer.com/.well-known/mcp-registry-auth; the private key lives outside the repo).
4. `mcp-publisher publish` from the repo root.
Bump `version` in `server.json` with each release.

## Other directories (free listings, do after the official one)
- Smithery (smithery.ai): submit the GitHub repo; stdio servers list fine.
- Glama (glama.ai/mcp/servers): auto-indexes GitHub repos with an MCP manifest; claim the listing.
- PulseMCP, mcp.so, Cursor directory: submit the repo URL.

## Claude Desktop extension (.mcpb)
Claude Desktop accepts one-click bundles. A `manifest.json` + the Python server packed with `npx @anthropic-ai/mcpb pack`
would let non-terminal users install it from the Desktop UI. Depends on bundling Python and the tools; do this
together with the desktop-app work.

## Privacy line for every listing
"Runs locally over stdio. Reads only the user's raw file on their machine. Returns derived results and single-site
lookups to the assistant; never the raw file. Health tables are separate tools fetched only on request."
