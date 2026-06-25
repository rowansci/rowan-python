# MSA

## Input

A list of protein sequences as amino-acid strings, passed as `initial_protein_sequences`. Pass one sequence for a single-chain MSA, or multiple sequences for a paired MSA. Generates and formats multiple-sequence-alignment data for downstream use with co-folding models.

## Example

```python
import rowan

folder = rowan.get_folder("examples")

wf = rowan.submit_msa_workflow(
    initial_protein_sequences=[
        "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
    ],
    output_formats={rowan.MSAFormat.BOLTZ},   # defaults to {rowan.MSAFormat.COLABFOLD}
    folder=folder,
)

result = wf.result()
result.download_files(rowan.MSAFormat.BOLTZ, path="msa_files")   # saves a .tar.gz archive
```

## Settings

- `output_formats` (default `{rowan.MSAFormat.COLABFOLD}`): which co-folding MSA formats to generate. Use `rowan.MSAFormat.BOLTZ` (for Boltz-1 / Boltz-2), `rowan.MSAFormat.CHAI` (for Chai-1), or `rowan.MSAFormat.COLABFOLD` (raw MMSeqs2 output for AlphaFold-derived pipelines). Pass a set.

## Retrieving output

The result does not expose alignment data as fields. Instead call `result.download_files(format, path=...)` to save the MSA files for a given format to disk. Passing a format that was not in `output_formats` raises `ValueError`. Omit `format` to download every requested format. Each download is a `.tar.gz` archive whose contents depend on the format: `boltz` gives per-sequence CSVs (`seq_0.csv`, `seq_1.csv`, ...), `chai` gives an `aligned.pqt` file, and `colabfold` gives `unpaired/` and `paired/` `.a3m` files.
