"""Example: submit a workflow with a webhook callback and manage webhook secrets."""

import rowan

# rowan.api_key = "rowan-sk..."

# Create a secret (idempotent — no-op if one already exists)
secret = rowan.create_webhook_secret()

# Rotate to a new secret (invalidates the old one)
# secret = rowan.rotate_webhook_secret()

oseltamivir = "C1CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H]([NH3+])C1CCC1"

workflow = rowan.submit_admet_workflow(
    initial_smiles=oseltamivir,
    name="Oseltamivir ADMET (webhook)",
    webhook_url="https://your-server.com/webhook",
)

print(f"Submitted: https://labs.rowansci.com/workflow/{workflow.uuid}")

# When the workflow completes, Rowan will POST the result to webhook_url with an
# X-Rowan-Signature header. Verify it on your server using HMAC-SHA256 with:
#
#   rowan.verify_webhook_secret(raw_body, signature_header, secret)
