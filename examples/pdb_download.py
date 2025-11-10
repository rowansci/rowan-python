from pathlib import Path

import rowan

# rowan.api_key = ""

proteins = rowan.list_proteins()
for protein in proteins:
    print(protein.name)
    protein.download_pdb_file(name=protein.name, path=Path("pdb_files"))
