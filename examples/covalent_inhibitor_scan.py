"""Submit a surface scan for prepared cathepsin K/S2 with matching atom references."""

import argparse
from pathlib import Path

import rowan


def main() -> None:
    """Upload a prepared complex and create a draft, or run it with `--run`."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pdb_path", type=Path, help="Prepared cathepsin K/S2 PDB with explicit hydrogens"
    )
    parser.add_argument("--two-dimensional", action="store_true", help="Include proton transfer")
    parser.add_argument("--run", action="store_true", help="Submit and wait for the calculation")
    args = parser.parse_args()

    cysteine = rowan.ResidueRef("polymer", "A", "A.25")
    ligand = rowan.ResidueRef("non_polymer", "A.300")
    coordinates = [
        rowan.ScanSettings(
            coordinate=rowan.DistanceCoordinate(
                atoms=(rowan.ResidueAtomRef(cysteine, "SG"), rowan.ResidueAtomRef(ligand, "C27"))
            ),
            start=1.6,
            stop=4.6,
            num=31,
            force_constant=100,
        )
    ]
    if args.two_dimensional:
        histidine = rowan.ResidueRef("polymer", "A", "A.162")
        hydrogen = rowan.ResidueAtomRef(histidine, "HD1")
        coordinates.append(
            rowan.ScanSettings(
                coordinate=rowan.LinearCombinationCoordinate(
                    coordinates=[
                        rowan.DistanceCoordinate(
                            atoms=(rowan.ResidueAtomRef(histidine, "ND1"), hydrogen)
                        ),
                        rowan.DistanceCoordinate(
                            atoms=(rowan.ResidueAtomRef(ligand, "N35"), hydrogen)
                        ),
                    ],
                    coefficients=[1, -1],
                ),
                start=-2,
                stop=2,
                num=21,
            )
        )
    region = [
        rowan.ResidueRef("polymer", "A", f"A.{number}")
        for number in (18, 19, 20, 23, 24, 25, 26, 27, 160, 161, 162, 163, 164, 182, 183, 184, 185)
    ]
    settings = rowan.CovalentInhibitorScanSettings(
        ligand=ligand, reactive_region=[*region, ligand], scan_coordinates=coordinates
    )
    protein = rowan.upload_protein("Prepared cathepsin K/S2", args.pdb_path)
    workflow = rowan.submit_covalent_inhibitor_scan_workflow(
        protein=protein,
        reactant_smiles="Cn1cnc2cnc(C#N)nc21",
        settings=settings,
        name="Cathepsin K/S2 covalent surface",
        is_draft=not args.run,
    )
    print(f"https://labs.rowansci.com/covalent-inhibitor-scan/{workflow.uuid}")
    if not args.run:
        return
    result = workflow.result()
    print(result)
    for window in result.windows:
        print(f"Grid index {window.index}: centers={window.centers}, samples={window.n_samples}")
    if result.surface is not None:
        print("Bin edges (Å):", result.surface.bin_edges)
        print("Free energies (kcal/mol):", result.surface.free_energies)
        print("Uncertainty (kcal/mol):", result.surface.uncertainty)
    for message in result.messages:
        print(f"{message.type}: {message.title}: {message.body}")


if __name__ == "__main__":
    main()
