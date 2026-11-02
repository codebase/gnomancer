# Development

```bash
./setup.sh            # venv, editable install, tools
.venv/bin/pytest
python -m build       # sdist + wheel; releases publish to PyPI from GitHub Actions (trusted publishing)
```

CI runs the tests and builds a demo report on every push. Real raw-data files and `*_analysis/` folders are
git-ignored; never commit them.

## Website

`site/` is the landing page at https://gnomancer.com, served as a Cloudflare Worker with static assets (`wrangler.toml`).
Deploy with `npx wrangler deploy` (Node 22+, `wrangler login` once). It is static: no cookies, no analytics.

## Release checklist

1. Bump the version in `pyproject.toml`, `gnomancer/__init__.py`, `gnomancer/mcp_server.py`, `server.json`, `.claude-plugin/plugin.json`.
2. `gh release create vX.Y.Z` → the publish workflow uploads to PyPI (trusted publishing).
3. `mcp-publisher login http --domain gnomancer.com --private-key <seed>` then `mcp-publisher publish` (key lives outside the repo).
4. Homebrew tap: update `url` and `sha256` in `codebase/homebrew-tap/Formula/gnomancer.rb` from the PyPI sdist (`https://pypi.org/pypi/gnomancer/X.Y.Z/json`).
5. If a pinned tool version changes, update `SHA256` in `gnomancer/setup_cmd.py`.
