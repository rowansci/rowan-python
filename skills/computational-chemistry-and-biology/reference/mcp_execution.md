# MCP execution

Authentication is attached to the MCP connection.

1. Start with `account_status`; use `mcp_supported_workflows` and available credits. Call
   `discover_workflow` for the exact schema—never guess parameters or enums.
2. For related real work, create or select a Rowan folder and pass its `folder_uuid` to each
   workflow. Import structures server-side with `import_structure`; inspect proteins, select one
   biologically relevant chain set rather than symmetry duplicates, and prepare apo receptors
   without heterogens unless a workflow needs an explicitly retained ligand or cofactor.
3. Create a draft first. Check its credit cap and dispatch estimate, then submit only within the
   user's authorization. Retain workflow UUIDs.
4. Retrieve results narrowly: use `wait_for_workflow_result` for known fields, otherwise inspect
   one preview and request selected fields with `result_path`, `offset`, and `limit`. Inspect or
   download one trajectory replicate before requesting all; fetch deferred files only when needed.
5. Stop invalid work before replacement; delete only terminal workflows when requested. Failed
   runs may charge credits. Preserve the UUID and use its diagnostics to choose a documented,
   scientifically meaningful retry rather than repeating the same failure.
6. If no named tool fits, search `discover_sdk_actions` narrowly and paginate when needed. Execute
   only advertised JSON-native actions; explain an unavailable file/object boundary rather than
   inventing an encoding.

Keep total spending within the authorized budget, avoid dumping large results into context, and
reuse discovered schemas for related runs.
