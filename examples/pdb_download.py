from pathlib import Path

import rowan

# Set your API key or use the ROWAN_API_KEY environment variable
# rowan.api_key = "rowan-sk..."

proteins = rowan.list_proteins()
for protein in proteins:
    print(protein.name)
    protein.download_pdb_file(name=protein.name, path=Path("pdb_files"))
