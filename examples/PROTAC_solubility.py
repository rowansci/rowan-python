# ruff: noqa: E501

import rowan

PROTACs = {
    4564: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCCN3C=C(CCCOC(=O)NCC4=CC=C(C(=O)NC5=CC=CC=C5N)C=C4)N=N3)C(C)(C)C)C=C2)SC=N1",
    2438: "CC(=O)N[C@H](C(=O)N1C[C@H](O)C[C@H]1C(=O)NCC1=CC=C(C2=C(C)N=CS2)C=C1OCCCCCCCCCCCC(=O)NC1=CC=C(C(=O)NC2=CC=CC=C2N)C=C1)C(C)(C)C",
    727: "NC1=CC=CC=C1NC(=O)C1=CC=C(NC(=O)CCCCCNC(=O)COC2=CC=CC3=C2C(=O)N(C2CCC(=O)NC2=O)C3=O)C=C1",
    728: "NC1=CC=CC=C1NC(=O)C1=CC=C(NC(=O)CCCCCCCCCCCNC(=O)COC2=CC=CC3=C2C(=O)N(C2CCC(=O)NC2=O)C3=O)C=C1",
    729: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    730: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    1843: "O=C(CCCCCCNC(=O)C1=CC=C(NC(=O)COCCOCCNC(=O)COC2=CC=CC3=C2C(=O)N(C2CCC(=O)NC2=O)C3=O)C=C1)NO",
    2110: "CC(C)C[C@H](NC(=O)[C@@H](O)[C@H](N)CC1=CC=CC=C1)C(=O)NCCCCCC(=O)NC1=CC=C(NC(=O)CCCCCCC(=O)NO)C=C1",
    2111: "CC(C)C[C@H](NC(=O)[C@@H](O)[C@H](N)CC1=CC=CC=C1)C(=O)NCCOCCOCCNC(=O)C1=CC=CC(C(=O)NO)=C1",
    2112: "CC(C)C[C@H](NC(=O)[C@@H](O)[C@H](N)CC1=CC=CC=C1)C(=O)NCCOCCOCC(=O)NC1=CC=CC(C(=O)NO)=C1",
    2419: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2420: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2421: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2422: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2423: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)COCCCCCCCCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2424: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)CCCCCCCCOCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2425: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)COCCCCCCOCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2426: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)COCCCCCCCCCOCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2427: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)COCCOCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
    2428: "CC1=C(C2=CC=C(CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@@H](NC(=O)COCCOCCOCC(=O)NC3=CC=C(C(=O)NC4=CC=CC=C4N)C=C3)C(C)(C)C)C=C2)SC=N1",
}

protac_solubility_folder = rowan.create_folder(name="PROTAC Solubility")

workflows = []
for id, smiles in PROTACs.items():
    workflows.append(
        rowan.submit_solubility_workflow(
            initial_smiles=smiles,
            solubility_method="fastsolv",
            solvents=["CS(=O)C"],
            temperatures=[293.15],
            folder_uuid=protac_solubility_folder.uuid,
            name=f"solubility {id}",
        )
    )


for workflow in workflows:
    workflow.wait_for_result()
    workflow.fetch_latest(in_place=True)

print(
    [
        (workflow.name, workflow.data["solubilities"]["CS(=O)C"]["solubilities"])
        for workflow in workflows
    ]
)
