from pathlib import Path

import rowan

# rowan.api_key = ""

pose_analysis_md_workflows = rowan.list_workflows(workflow_type="pose_analysis_md", status=2)
for workflow in pose_analysis_md_workflows:
    print(workflow.name)
    workflow.download_dcd_files(replicates=[0], name=workflow.name, path=Path("dcd_files"))
