from pathlib import Path

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

pose_analysis_md_workflows = rowan.list_workflows(workflow_type="pose_analysis_md", status=2)
for workflow in pose_analysis_md_workflows:
    print(workflow.name)
    workflow.download_dcd_files(replicates=[0], name=workflow.name, path=Path("dcd_files"))
