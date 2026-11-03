import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _fake_report(d: Path):
    d.mkdir()
    (d / "summary.json").write_text(json.dumps(dict(input=str(d / "nofile.txt"), args=[], glance={"who": ["Y: R-M269"]}, source="AncestryDNA",
                                                   count=10, qc={"sex_inferred": "Male"}, y={"haplogroup": "R-M269"}, mt={"haplogroup": "H1"},
                                                   combos=[["ABO blood group", "O", ""]], clinvar={"hits": []}, pharmcat={"calls": []})))
    (d / "REPORT.md").write_text("# Report\n\n## 1. Input\n\nvendor\n\n## 4. Paternal line (Y chromosome)\n\n**Haplogroup R-M269**\n\n## 8. Marker panel\n\ntable\n")
    (d / "marker-panel.tsv").write_text("rsid\tcategory\tgene\tlabel\tref\teffect\tchrom\tpos\tgenotype\teffect_count\tinterpretation\tnote\n"
                                        "rs762551\ttrait\tCYP1A2\tCaffeine\tA\tC\t15\t1\tAA\t0/2\tfast metabolizer\t\n")
    (d / "summary-prompt.md").write_text("PROMPT")


def test_tools_read_a_report(tmp_path):
    from gnomancer import mcp_server as m
    d = tmp_path / "rep"
    _fake_report(d)
    s = m.summary(str(d))
    assert s["y"]["haplogroup"] == "R-M269" and s["at_a_glance"]["who"] == ["Y: R-M269"]
    assert "Paternal line" in m.section(None, str(d))
    assert "R-M269" in m.section("paternal", str(d))
    with pytest.raises(ValueError):
        m.section("nonexistent", str(d))
    r = m.lookup("rs762551", str(d))
    assert r["genotype"] == "AA" and r["interpretation"] == "fast metabolizer"
    r = m.lookup("rs999999", str(d))
    assert r["genotype"] is None
    assert m.narrative_prompt(str(d)) == "PROMPT"
    assert m.clinvar_findings(str(d)) == {"hits": []}
    with pytest.raises(ValueError):
        m.summary(str(tmp_path / "missing"))


def test_stdio_server_lists_tools():
    """Spawn the real server over stdio and list tools with the MCP client SDK."""
    import anyio
    from mcp import ClientSession
    from mcp.client.stdio import stdio_client, StdioServerParameters

    async def go():
        params = StdioServerParameters(command=sys.executable, args=["-m", "gnomancer", "serve"], cwd=str(ROOT))
        async with stdio_client(params) as (r, w):
            async with ClientSession(r, w) as s:
                await s.initialize()
                tools = {t.name for t in (await s.list_tools()).tools}
                res = {str(x.uri) for x in (await s.list_resources()).resources}
                return tools, res

    tools, res = anyio.run(go)
    assert {"analyze", "summary", "section", "lookup", "clinvar_findings", "pharmacogenomics", "narrative_prompt", "set_narrative"} <= tools
    assert "gnomancer://guide/get-your-data" in res
