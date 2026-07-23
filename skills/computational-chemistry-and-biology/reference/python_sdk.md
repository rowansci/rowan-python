# Python SDK

Use this path when working through the Rowan Python SDK. Use the resulting environment's Python for every later step.

## Setup

**1. Install `rowan`** (PyPI: `rowan-python`). Skip if `python3 -c "import rowan"` already succeeds. Otherwise add it through the project's environment manager so it survives syncs. A bare `pip install` into a managed environment gets wiped:

- **uv** (`uv.lock`): `uv add rowan-python`, then run with `uv run python`
- **pixi** (`pixi.toml`): `pixi add --pypi rowan-python`, then `pixi run python`
- **poetry** (`poetry.lock`): `poetry add rowan-python`, then `poetry run python`
- **plain venv or system**: `python3 -m pip install rowan-python`

**2. Set the key and verify.** Put `ROWAN_API_KEY=your-key-here` in a project `.env` (gitignored, and create a key at <https://labs.rowansci.com/account>), then run the preflight. It finds the nearest `.env`, validates the key against the API, and prints remaining credits. Run it from the project root, where `.env` lives, and point it at this skill's base directory rather than changing into the skill:

`python3 "<skill_dir>/scripts/check_env.py"`

It exits with actionable guidance for a missing key, a rejected or expired key, or an unreachable API. Since shell state does not persist between commands, load `.env` in the same command as every later call, such as `[ -f .env ] && set -a && source .env && set +a; <python ...>`, or set `rowan.api_key = os.environ["ROWAN_API_KEY"]` in Python.

## Run and retrieve

Every workflow is submitted and retrieved the same way. `submit_*_workflow(...)` returns a `Workflow` immediately, and results are fetched separately.

```python
import rowan

folder = rowan.get_folder("my-project")
wf = (
    rowan.submit_
    < workflow
    > _workflow(
        ...,  # use the scientific arguments selected by the entry skill
        folder=folder,
        name="my run",
    )
)
result = wf.result()
```

**Universal submit arguments.** Every `submit_*_workflow` accepts these in addition to its scientific arguments:

- `folder=` or `folder_uuid=` sets where the run lands. `rowan.get_folder("a/b/c")` returns the folder, creating the path if missing. For a specific project or benchmark, make a dedicated folder and route all its runs there.
- `name=` sets the label shown in the Rowan UI.
- `max_credits=` is an optional runtime cap. Omit it unless the user requests one or a bounded budget or unusual cost risk calls for it. A low cap can permanently stop a run.
- `webhook_url=` makes Rowan POST the result there on completion.
- `is_draft=True` stages the run without starting it.

**Retrieve results:**

- **Block:** `result = wf.result()` waits and returns a typed result. It raises `rowan.WorkflowError` if the run failed or stopped.
- **Partial or progress:** `wf.result(wait=False)` returns what is ready now. `for result in wf.stream_result():` yields partial results and then the final result.
- **Fire-and-forget:** preserve `wf.uuid`, then reconnect later with `rowan.retrieve_workflow(uuid).result()`. `rowan.list_workflows()` lists recent runs. The UUID from `labs.rowansci.com/calculation/<uuid>` is a workflow UUID; use `retrieve_workflow`, not `retrieve_calculation`.

Use inline submit-and-wait for quick questions or one result. For long-running or multi-step experiments, write a reproducible script and prefer fire-and-forget retrieval or a webhook over blocking the session.

**Estimate before committing:**

```python
draft = rowan.submit_ < workflow > _workflow(..., folder=folder, is_draft=True)
print(draft.dispatch_info())
draft.submit_draft()  # or draft.delete()
```

**Status and control:** `wf.done()` and `wf.get_status()` check state without waiting, `wf.stop()` cancels a running workflow, and `wf.delete()` removes it and its data. `wf.update(...)` edits workflow metadata (`name`, folder, `starred`, `public`, `email_when_complete=True` to email on completion).

**Resubmit:** call `rowan.submit_workflow(workflow_type, workflow_data, initial_molecule=..., folder=folder)`. Pass a prior run's `result.data` unchanged for an identical rerun, edit it to change settings, or use a matching workflow type and data to submit another workflow.

**Perturb geometry before resubmission:** `mol.perturb()` adds small Gaussian noise. `mol.displace_along_mode(mode, displacement)` shifts atoms along a vibrational mode; for a transition state, use the imaginary mode from a calculation that included frequencies.
