# Rowan Python Library

[![pypi](https://img.shields.io/pypi/v/rowan-python.svg)](https://pypi.python.org/pypi/rowan-python)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://docs.astral.sh/uv/)
[![ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)](https://github.com/charliermarsh/ruff)
[![Downloads](https://img.shields.io/pypi/dm/rowan-python.svg)](https://pypi.python.org/pypi/rowan-python/)
[![License](https://img.shields.io/github/license/rowansci/rowan-python)](LICENSE)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/rowansci/rowan-python/test.yml?branch=master&logo=github-actions)](https://github.com/rowansci/rowan-python/actions)
<!-- Enable these badges with the corresponding tooling/services.
[![Markdown style: rumdl](https://img.shields.io/badge/md%20style-rumdl-000000.svg)](https://rumdl.dev)
[![Typing: ty](https://img.shields.io/badge/typing-ty-EFC621.svg)](https://github.com/astral-sh/ty)
[![Codecov](https://img.shields.io/codecov/c/github/rowansci/rowan-python)](https://codecov.io/gh/rowansci/rowan-python)
-->



The Rowan Python library provides convenient access to the Rowan API from applications written in the Python language.

## Documentation

The documentation is available [here](https://docs.rowansci.com/python-api).

## Agent skill

Ships with a [computational chemistry and biology skill](skills/computational-chemistry-and-biology/)
that helps coding agents choose and run Rowan workflows through either Rowan MCP tools or the
Rowan Python SDK.

### Claude Code

```bash
claude plugin marketplace add https://github.com/rowansci/rowan-python.git
claude plugin install computational-chemistry-and-biology@rowan
```

### Codex

```bash
codex plugin marketplace add rowansci/rowan-python --ref master
codex plugin add computational-chemistry-and-biology@rowan
```

### Rowan MCP server (optional)

The skill also works with the hosted Rowan MCP server, which provides Rowan tools directly to the
agent. Add it once per client:

```bash
# Claude Code
claude mcp add --transport http rowan https://mcp.rowansci.com/

# Codex
codex mcp add rowan --url https://mcp.rowansci.com/
codex mcp login rowan
```

Both clients prompt for OAuth sign-in on first use. Without the server, the skill falls back to the
Rowan Python SDK, which authenticates with `ROWAN_API_KEY` instead.

Start a new Claude Code or Codex session after installation. For manual installation, download the
[latest skill ZIP](https://github.com/rowansci/rowan-python/releases/download/skill-latest/computational-chemistry-and-biology-skill.zip)
and extract it into your agent's skills directory.

## Running examples

To run the examples, you will need to set your ROWAN_API_KEY environment variable or set it directly in the script.
If running the examples in a cloned version of the repository, you can add your api key to a local `.env` file, which will automatically be loaded into the environment by direnv (if installed).


## Issues

To report issues, please use the "Issues" tab above.

*Corin Wagen, 2023*
