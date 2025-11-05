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
    output_formats=[MSAFormat.BOLTZ],
    name="Boltz Paired MSA Example",
)

msa_workflow.wait_for_result().fetch_latest(in_place=True)

msa_workflow.download_msa_files(MSAFormat.BOLTZ, path=msa_directory)

tar_path = next(msa_directory.glob("*.tar.gz"))
with tarfile.open(tar_path, "r") as tar_ref:
    tar_ref.extractall(msa_directory)

tar_path.unlink()

# by default, the csv files are named seq_<index>.csv
# to use them for boltz, you need to set the msa path in the protein section of the boltz input
# yaml to the path of the relevant csv file
