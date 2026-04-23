import tarfile
from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."
folder = rowan.get_folder("examples")

msa_directory = Path("msa_directory")

msa_workflow = rowan.submit_msa_workflow(
    initial_protein_sequences=[
        "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
    ],
    output_formats=["boltz"],
    name="Boltz MSA Example",
    folder=folder,
)

print(f"View workflow privately at: https://labs.rowansci.com/msa/{msa_workflow.uuid}")

msa_result = msa_workflow.result()

msa_result.download_files("boltz", path=msa_directory)

tar_path = next(msa_directory.glob("*.tar.gz"))
with tarfile.open(tar_path, "r") as tar_ref:
    tar_ref.extractall(msa_directory)

tar_path.unlink()

# by default, the csv files are named seq_<index>.csv
# to use them for boltz, you need to set the msa path in the protein section of the boltz input
# yaml to the path of the relevant csv file
