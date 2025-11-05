import tarfile
from pathlib import Path

# from chai_lab.chai1 import run_inference
from stjames import MSAFormat

import rowan

# rowan.api_key = ""
example_fasta = (
    ">protein|name=example-protein\n"
    "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW\n"
)
fasta_path = Path("/tmp/input.fasta")
fasta_path.write_text(example_fasta)

output_dir = Path("/tmp/outputs_2")
msa_directory = Path("msa_directory")

msa_workflow = rowan.submit_msa_workflow(
    initial_protein_sequences=[
        "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
    ],
    output_formats=[MSAFormat.CHAI],
    name="CHAI MSA Example",
)

msa_workflow.wait_for_result().fetch_latest(in_place=True)

msa_workflow.download_msa_files(MSAFormat.CHAI, path=msa_directory)

tar_path = next(msa_directory.glob("*.tar.gz"))
with tarfile.open(tar_path, "r") as tar_ref:
    tar_ref.extractall(msa_directory)

tar_path.unlink()


# run_inference(fasta_file=fasta_path,
#         output_dir=output_dir,
#         num_trunk_recycles=3,
#         num_diffn_timesteps=200,
#         seed=42,
#         device="cpu", # or "cuda:0"
#         use_esm_embeddings=True,
#         use_msa_server=False,
#         use_templates_server=False,
#         msa_directory=msa_directory)
