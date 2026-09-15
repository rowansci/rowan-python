# AI Agent Development Guide

This document provides essential guidance for AI agents working on this repository.

## AI skills

Shared development conventions live in `.agents/skills/`. Read `write-code`,
`write-docstrings`, and `write-tests` when working on the corresponding code.
These files and this guide are tracked development resources, excluded from Python distributions.
Keep personal guidance in global agent settings or locally excluded files.

The cookiecutter migration is staged: use Google-style docstrings and keep existing lint exceptions;
use mypy until the ty migration. Markdown checks and Codecov uploads remain disabled.
Do not apply the deferred conventions as a repository-wide cleanup.

## Repository overview

rowan-python is the Python SDK for the Rowan computational chemistry platform. It wraps the stjames data model library and provides user-facing workflow submission and result retrieval.

Structure:
- `rowan/` - source code (flat layout)
- `rowan/workflows/` - workflow submission functions and result types
- `examples/` - usage examples for each workflow
- `docs/` - documentation (mkdocs)
- `.github/workflows/` - CI/CD configuration

Key relationships:
- **stjames** (`../stjames`) - data model dependency, imported as `stjames`. Defines workflow models, settings, validation, and engine compatibility.
- **rowan-python** - public version of this repo. Push via `git push public master` from this repo.

Python version: >=3.12

## Essential commands

```bash
# Setup
uv sync --locked               # Install dependencies
uv run prek install --hook-type pre-commit --hook-type pre-push

# Code quality (these are the pre-commit hooks)
uv run ruff format .           # Format code
uv run ruff check . --fix      # Lint code
uv run mypy .                  # Type check

# Testing
uv run pytest                  # Run tests and doctests
uv run prek run --all-files    # Run pre-commit checks

# Documentation
uv run --group docs mkdocs serve
uv run --group docs mkdocs build --clean

# Run a specific example
uv run python examples/basic_calculation.py
```

## Before every commit

- Run `uv run ruff format .`, `uv run ruff check . --fix`, `uv run mypy .`
- Pre-commit hooks (`prek.toml`) run these automatically on commit
- Tests run on pre-push and in CI

## Code conventions

### Documentation updates

- Treat requests to check or update docs as keeping existing descriptions, examples, docstrings, and skill references accurate.
- Make corrections where the relevant behavior is already documented. Don't add feature announcements, introductory callouts, or new sections just because a feature changed.
- If a document doesn't discuss the affected behavior, usually leave it alone. Preserve its existing scope and emphasis unless the user asks for expanded documentation.

### Docstrings

Follow the [docstring skill](.agents/skills/write-docstrings/SKILL.md) for Google-style
formatting and description conventions.

### Type annotations

- All functions must have complete type annotations
- Modern syntax: `list[str]`, `dict[str, int]`, `str | None`
- Import from `typing` only when necessary

### Code formatting

Via ruff:
- Line length: 100
- Indentation: 4 spaces

### Imports

- Absolute imports preferred
- No wildcard imports except in `__init__.py`
- Import sorting handled by ruff (isort)

## Workflow development guidelines

### Validation
- Use stjames validation when possible. Don't duplicate validation that stjames model validators already handle (e.g. engine/method compatibility, solvent checks). Only add rowan-side validation when stjames doesn't cover it.

### stjames type aliasing
- Users should never need to `import stjames` directly. All user-facing stjames types must be aliased in `rowan/__init__.py`. If an example requires a stjames import, that's a signal the type needs to be aliased.

### Serialization
- Use `serialize_as_any=True` on `model_dump` when the workflow has union-typed fields (e.g. `ConformerGenSettingsUnion`, `MultiStageOptSettings` containing `Settings` subfields). Without it, pydantic may silently drop subclass-specific fields during serialization.

### Defaults
- When hardcoding a default value for a stjames field, make sure it matches the corresponding default in the stjames workflow model.

### Type hints in function signatures
- Don't use pydantic-specific types (`PositiveInt`, `NonNegativeInt`, etc.) in plain function signatures. They don't enforce constraints outside pydantic models and are misleading. Use plain `int`, `float`, etc. — stjames validates downstream when the model is constructed.

### Testing
- After editing a workflow, run its relevant example before committing (not on every change) to catch breakage early.

## Git authorization policy

- Prepare changes on a local feature branch. Local commits are allowed as part of authorized development work.
- The user handles pushing and opening PRs. Do not push branches, create PRs, or change remote branches unless the user explicitly requests that specific remote action.
- A request to get changes ready for a PR means prepare the local branch, not publish it or open a PR.

**Never add yourself as a commit author or co-author.** Do not include `Co-Authored-By:`, `Author:`, or any similar trailer that attributes the commit to an AI model or tool. Commits are attributed solely to the human developer.

## CI/CD

File: `.github/workflows/test.yml`

Triggers: all PRs, pushes to `master`

Checks:
1. `uv run ruff format --check --diff .` - format check
2. `uv run ruff check .` - lint check
3. `uv run mypy .` - type check
4. `uv run pytest --cov --cov-report=xml` - tests, doctests, and coverage
5. `uv run prek run --all-files check-toml check-yaml trailing-whitespace end-of-file-fixer` - file checks

Matrix: Python 3.12 and 3.14, ubuntu-latest

## Additional resources

- uv documentation: https://docs.astral.sh/uv/
- ruff documentation: https://docs.astral.sh/ruff
- pytest documentation: https://docs.pytest.org
