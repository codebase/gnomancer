# Repository rules for AI assistants

- Never add attribution trailers to commits, PRs, or release notes: no `Co-Authored-By: Claude`, no "Generated with Claude Code", nothing of the kind. A `commit-msg` hook in `.githooks/` enforces this; run `git config core.hooksPath .githooks` after cloning.
- Never commit raw DNA files, reports, or anything under `*_analysis/`.
- Health wording is hedged and sourced; see CONTRIBUTING.md.
