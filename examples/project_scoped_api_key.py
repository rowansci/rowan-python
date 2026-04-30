"""Example: mint per-project API keys for several projects at once.

A project-scoped key can only access the project it was created for —
useful for handing limited access to collaborators, CI jobs, or analyses
without exposing your full account.

The plaintext key is only returned once at creation time — store it
immediately. The caller must currently authenticate with an unscoped key
that has ``can_manage_api_keys`` permission.
"""

import rowan

# rowan.api_key = "rowan-sk..."

campaigns = ["CDK2 campaign", "BTK campaign", "EGFR campaign"]

results = []
for campaign in campaigns:
    project = rowan.create_project(name=campaign)
    created = rowan.create_api_key(
        name=f"{campaign.lower().replace(' ', '-')}-key",
        scope="read_write",
        valid_days=90,
        scoped_project_uuid=project.uuid,
    )
    results.append((project.name, created.api_key.name, created.key))

print("Plaintext keys (save these now — only returned once):")
for project_name, key_name, key in results:
    print(f"  {project_name} | {key_name} | {key}")
