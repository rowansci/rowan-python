---
name: cleanup-code
description: Polishes recently written or modified Python code for style, dead code, docstrings, type annotations, and checks. Use when coding task is completed.
allowed-tools:
  - Read
  - Glob
  - Grep
  - Edit
  - Bash
---

# Clean Up Code

Polish only the code changed for the current task, within $ARGUMENTS when specified.
Do not expand a tooling migration into a package-wide style, docstring, or test rewrite.
Follow the staged migration in [AGENTS.md](../../../AGENTS.md); deferred checks are
not requirements for this cleanup.

Use the `write-code` skill to understand code conventions.

1. Re-read every file that was created or modified.
2. Check for Pythonic style:
   - Use comprehensions instead of loops where clearer
   - Use f-strings instead of `.format()` or `%`
   - Use `pathlib` instead of `os.path`
   - Use dataclasses or named tuples for structured data
3. Remove dead code, unused imports, and obsolete commented-out code; preserve the intentionally deferred migration settings
4. Remove temporary tests introduced during the current task that provide no lasting value.
5. Verify docstrings on changed public APIs (follow `write-docstrings` skill conventions).
6. Verify complete type annotations on changed public functions and methods.
7. Run the full check suite:
   - `ruff format .`
   - `ruff check .`
   - `uv run mypy .`
   - `pytest`
8. All checks must pass before finishing.
