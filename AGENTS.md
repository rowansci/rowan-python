# AI Agent Development Guide

This document provides essential guidance for AI agents working on this repository.

## Repository overview

rowan-python-internal is the private Python SDK for the Rowan computational chemistry platform. It wraps the stjames data model library and provides user-facing workflow submission and result retrieval.

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
pixi install                    # Install dependencies

# Code quality (these are the pre-commit hooks)
pixi run fmt                    # Format code (ruff format)
pixi run lint                   # Lint code (ruff check --fix)
pixi run types                  # Type check (mypy)

# Testing
pixi run test                   # Run tests (pytest with doctests)
pixi run all                    # Run fmt + lint + types + test

# Run a specific example
pixi run python examples/basic_calculation.py
```

## Before every commit

- Run `pixi run fmt`, `pixi run lint`, `pixi run types`
- Pre-commit hooks (`.pre-commit-config.yaml`) run these automatically on commit
- No pytest in pre-commit hooks — tests run in CI

## Code conventions

### Docstrings

Format: reStructuredText-style. No types in docstrings, no leading articles.

```python
def process_data(input_data: list[str], threshold: int = 10) -> dict[str, int]:
    """Process input data and return summary statistics.

    :param input_data: strings to process
    :param threshold: minimum count threshold for inclusion
    :returns: mapping of categories to counts
    """
```

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

**Do not run `git commit`, `git push`, or any other git write commands without explicit permission from the user.**

**Never add yourself as a commit author or co-author.** Do not include `Co-Authored-By:`, `Author:`, or any similar trailer that attributes the commit to an AI model or tool. Commits are attributed solely to the human developer.

## CI/CD

File: `.github/workflows/test.yml`

Triggers: all PRs, pushes to `master`

Checks:
1. `pixi run fmt` - format check
2. `pixi run lint` - lint check
3. `pixi run types` - type check

Matrix: Python 3.14, ubuntu-latest

## Additional resources

- pixi documentation: https://pixi.sh
- ruff documentation: https://docs.astral.sh/ruff
- pytest documentation: https://docs.pytest.org
