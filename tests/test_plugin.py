"""Test the installable Rowan skill plugin contract."""

import json
import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).parents[1]
SKILL_ROOT = REPOSITORY_ROOT / "skills" / "computational-chemistry-and-biology"


def test_plugin_manifests_have_consistent_identity_and_version() -> None:
    """Keep Claude, Codex, and marketplace plugin identities consistent."""
    claude = json.loads((SKILL_ROOT / ".claude-plugin" / "plugin.json").read_text())
    codex = json.loads((REPOSITORY_ROOT / ".codex-plugin" / "plugin.json").read_text())
    marketplace = json.loads((REPOSITORY_ROOT / ".claude-plugin" / "marketplace.json").read_text())

    assert claude["name"] == codex["name"] == marketplace["plugins"][0]["name"]
    assert claude["version"] == codex["version"]
    assert marketplace["plugins"][0]["strict"] is True
    assert claude["skills"] == "./"


def test_skill_markdown_links_stay_inside_the_installable_plugin() -> None:
    """Keep skill references available after hosts copy the plugin directory."""
    markdown_files = [SKILL_ROOT / "SKILL.md", *sorted((SKILL_ROOT / "reference").glob("*.md"))]
    for markdown_file in markdown_files:
        for target in re.findall(r"\[[^]]+\]\(([^)]+\.md)\)", markdown_file.read_text()):
            resolved = (markdown_file.parent / target).resolve()
            assert resolved.is_relative_to(SKILL_ROOT.resolve())
            assert resolved.is_file(), f"{markdown_file}: missing {target}"
