"""Example: temporarily share a workflow with anyone who has its link."""

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable.
# rowan.api_key = "rowan-sk..."

workflow = rowan.retrieve_workflow("your-workflow-uuid")

# Temporary shares can last for up to 120 minutes. This does not make the workflow
# permanently public: `workflow.public` remains unchanged.
workflow = workflow.temporarily_share(duration_minutes=60)

print(f"Share URL: https://labs.rowansci.com/workflow/{workflow.uuid}")
print(f"Public until: {workflow.public_until}")
print(f"Temporary share active: {workflow.is_temporarily_public}")

# End access before the expiration time if needed.
# workflow = workflow.end_temporary_share()
