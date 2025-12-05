from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path


import ase.io
import ase.neighborlist


PATH_XYZ_CYCLOHEXANE = Path(__file__).parent / "files" / "cyclohexane.xyz"
PATH_XYZ_AZULENE = Path(__file__).parent / "files" / "azulene.xyz"
PATHS_XYZ = [PATH_XYZ_CYCLOHEXANE, PATH_XYZ_AZULENE]


# def test_inference_bonds_rdkit():
#     """
#         Tries to infer bonds from an .xyz using rdkit.
#     """
atoms = ase.io.read(PATH_XYZ_AZULENE, index=-1)


cutoffs = ase.neighborlist.natural_cutoffs(atoms)
nl = ase.neighborlist.NeighborList(
    cutoffs=cutoffs,
)
nl.update(atoms)


for path in PATHS_XYZ:
    chemgraph = cg.from_file(path_or_file=path, fmt="xyz")
    chemgraph = chemgraph.infer_bonds(method="rdkit")

    a = 1
