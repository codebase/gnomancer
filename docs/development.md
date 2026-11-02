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
