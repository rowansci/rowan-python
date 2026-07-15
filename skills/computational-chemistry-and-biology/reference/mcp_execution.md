# MCP execution

Use this path when Rowan MCP tools are available. Authentication is already attached to the MCP connection.

1. Call `account_status` before reading workflow references. Use `mcp_supported_workflows` to limit selection and report available user and organization credits.
2. Apply the scientific workflow guidance already selected by the entry skill. Retain its parameter rationale, interpretation guidance, caveats, and useful Python examples, but use `discover_workflow` for the exact MCP parameter contract.
3. Call `discover_workflow` with the enabled workflow slug. Treat the returned parameter specification as authoritative; never guess an exact field name, type, default, or enum.
4. Build a JSON `parameters` object from the user's inputs and the discovered specification. Ask only for scientifically meaningful missing values. For an attached PDB, SDF, MOL, MOL2, XYZ, EXTXYZ, or CIF file, call `import_structure` and use its `workflow_input` or returned molecule data in the discovered parameter. Inline import is for small text files; do not place a large file into model context.
5. Call `create_workflow_draft`. This validates the inputs and creates a non-running draft. Inspect its `dispatch_estimate`, `max_credits`, and `available_credits`. If the user has not granted standing authorization that covers this workflow, its parameters, and its per-run and cumulative credit cost, review the draft with them and request approval.
6. After specific draft approval, or when the draft remains within the user's standing authorization, call `submit_workflow_draft` with its `workflow_uuid`. Stop and request approval if its scope or cost exceeds that authorization. Preserve the UUID so the workflow can be recovered in another session.
7. Call `workflow_status` without waiting, or set `wait_seconds` to at most 60 for bounded polling. Do not repeatedly occupy a turn waiting for a long workflow. Failed or stopped workflows include a bounded `log_tail`; report the actionable failure without dumping internal logs.
8. Call `get_workflow_result` with no fields first. Inspect `available_fields` and `result_preview`, then request at most five desired fields per call. For nested arrays or objects, use `result_path`, `offset`, and `limit` rather than requesting the complete value.
9. Use `get_structure_file` for PDB, XYZ, SDF, MOL, or MOL2 output from a protein or calculation reference. Use `get_workflow_file` for MSA and trajectory archives. These return deferred resources so file content is loaded only when the user or host needs it.

Pass JSON-native parameters through MCP and use `import_structure` for supported structure files. If discovery requires another Rowan object and does not describe a supported JSON representation or reference, explain that limitation rather than inventing an encoding. A workflow reference remains useful when its executable example is Python: translate its scientific choices through the discovered MCP contract.
