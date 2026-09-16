---
name: write-code
description: Write or review code in this repo. Use when writing new code, reviewing existing code, or editing code.
---

# Writing Code

Use this skill when writing, reviewing, or editing code.

## Current migration stage

Follow the staged migration in [AGENTS.md](../../../AGENTS.md). Apply style guidance
within the requested work; do not rewrite existing code to satisfy deferred rules.
Keep existing lint exceptions and commented-out migration settings. Use ty, Rumdl, and the enabled
Ruff rules. Use Google-style docstrings.

## Purpose

The goal of this document is to provide guidance when coding in this repository. All code should be pythonic and easy to read.

Python version: see `requires-python` in `pyproject.toml`

## Before every commit

- Ensure all code has type annotations
- Use Google-style docstrings (NO types, NO leading articles)
- Run checks: `prek -a`
- Prek hooks will run automatically and must pass

## Code conventions

### Docstrings

See skill `write-docstrings`

### Type annotations

- All functions must have complete type annotations
- Import types from `typing` only when necessary (prefer built-ins)
- Use modern syntax: `list[str]`, `dict[str, int]` (not `List[str]`, `Dict[str, int]`)
- Do not use a bare `dict`, always annotate the type of the `dict` (e.g. `dict[str, float]`)
- Use `Any` in `dict` annotations only if absolutely necessary
- Use `|` for union types: `str | None`
- Avoid `from __future__ import annotations`
- Use modern numpy type hints where useful, and use TypeAlias to make code more readable, e.g.
  - `type Matrix[T: np.generic] = np.ndarray[tuple[int, int], np.dtype[T]]`

### Code formatting

Via ruff

- Line length: 100
- Indentation: 4 spaces (no tabs except Makefiles)

### Naming conventions

- Functions/methods: snake_case
- Variables: snake_case
- Constants: UPPER_SNAKE_CASE
- Classes: PascalCase
- Modules: snake_case
- Private attributes/methods: _leading_underscore

### Imports

- Absolute imports preferred
- Group imports: standard library, third-party, local
- No wildcard imports (`from module import *`) except in `__init__.py`
- Import sorting handled by ruff (isort)
- Do not use import statements within a function unless absolutely necessary

### Error handling

- Prefer specific exceptions; define custom exceptions for domain errors
  (e.g. `SpamTypeError(ValueError)`, `SpamEatingError(RuntimeError)`)
- Avoid bare `except:` clauses

### General style

- Use f-strings for string formatting
- Prefer list/dict comprehensions over loops when appropriate
- Use `pathlib.Path` for file operations instead of `os.path`
- Use dataclasses and prefer the settings `slots=True` and `frozen=True`
- Don't use `.0` to indicate floats
- Prefer `strict=True` in `zip` and `itertools.batched`
- Use a guard case in all `match`/`case` statements (i.e. `case _:`)
- Don't use the filename `types.py` as it conflicts with a builtin library
- Avoid using `cast` to achieve the correct type

## Essential commands

Code quality tools (ruff, ty, pytest) are configured per-package in `pyproject.toml`.

```bash
# Setup
prek install                    # Install git hooks

# Code quality
ruff format .                   # Format code
ruff check .                    # Lint code
uv run ty check                 # Type check
prek -a                         # Run all prek hooks
prek run <hook-id>              # Run specific hook

# Testing
pytest                          # Run tests
pytest --cov                    # Run tests with coverage

# Package management
uv sync                         # Install dependencies
uv add <package>                # Add dependency
uv add --dev <package>          # Add dev dependency
```

## Direnv integration

Direnv auto-activates virtual environments when entering package directories.
Run `direnv allow .` if you see permission errors.

## Testing

- Use pytest; tests live in `tests/`; doctests are auto-discovered in source
- Use standard pytest format (e.g. don't use classes to hold tests)
- For float comparisons, use `approx` or `assert_almost_equal` (imported as `aae`), prefer default thresholds

See the skill `write-tests` for more detail, but only if actively writing tests

## CI/CD

CI runs non-mutating Rumdl and Ruff checks, ty, file checks, and pytest with coverage. Run `prek -a`
to reproduce most of the set locally.

## Git development guidelines

When working on a new feature, bugfix, etc.
If Claude is already working in this directory, offer to start a new worktree.
Create a new branch for all new work (avoid developing on master)

e.g.
`git worktree add ../{package_name}-{feature_name} && cd ../{package_name}-{feature_name}` # optional
`git switch -c feat/jalapeno_spam`

Break up work into manageable commits, keeping track of what changed in each commit.
Prepare changes on a local feature branch. Do not push or open a PR unless the user
explicitly requests that remote action.

If a worktree was created, offer to clean it up.

### Commit guidelines

Format: conventional commits recommended

- `feat:` - new features
- `fix:` - bug fixes
- `docs:` - documentation changes
- `test:` - test changes
- `refactor:` - code refactoring
- `chore:` - maintenance tasks

Never attribute commits or co-author trailers to an AI model or tool.

#### Example

```bash
git commit -m "feat: add Jalapeno Spam

- Adds Jalapeno Spam
- Tests that it is spicy"
```

### Documentation

Use skill `write-docstrings`

## Additional resources

- [uv](https://docs.astral.sh/uv)
- [ruff](https://docs.astral.sh/ruff)
- [ty](https://github.com/astral-sh/ty)
- [pytest](https://docs.pytest.org)
- [prek](https://prek.j178.dev)
- [Google docstring style](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings)

## Troubleshooting

### Prek hook failures

Formatting issues:

- Usually auto-fixed by ruff
- Re-stage files: `git add .`
- Try committing again

Lint issues:

- Read error message for specific rule
- Fix; if `# noqa: <rule>` if absolutely necessary, ask before adding

Type issues:

- Add missing type annotations
- Fix type mismatches
- Use `uv run ty check` to verify locally

Test failures:

- Fix failing tests or code
- Run `pytest -v` for detailed output
- Run specific test: `pytest tests/test_file.py::test_name`

## Miscellaneous

Do not add new packages without explicitly asking first.
Never add ignores to the pyproject.toml formatting, linting, or type checking without explicitly asking
