import tarfile
from pathlib import Path

from stjames import MSAFormat

import rowan

# rowan.api_key = ""

msa_directory = Path("msa_directory")

msa_workflow = rowan.submit_msa_workflow(
    initial_protein_sequences=[
        "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
        "VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
    ],
    output_formats=[MSAFormat.COLABFOLD],
    name="Colabfold Paired MSA Example",
)

msa_workflow.wait_for_result().fetch_latest(in_place=True)

msa_workflow.download_msa_files(MSAFormat.COLABFOLD, path=msa_directory)

tar_path = next(msa_directory.glob("*.tar.gz"))
with tarfile.open(tar_path, "r") as tar_ref:
    tar_ref.extractall(msa_directory)

tar_path.unlink()

# This produces two folders, one with unparied msas called unpaired and one with paired msas
# called paired
