from pathlib import Path

import rowan

# Set ROWAN_API_KEY environment variable to your API key or set rowan.api_key directly
# rowan.api_key = "rowan-sk..."

proteins = rowan.list_proteins()
for protein in proteins:
    print(protein.name)
    protein.download_pdb_file(name=protein.name, path=Path("pdb_files"))
