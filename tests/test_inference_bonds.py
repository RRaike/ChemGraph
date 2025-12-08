from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path

PATH_XYZ_CYCLOHEXANE = Path(__file__).parent / "files" / "cyclohexane.xyz"
PATH_XYZ_AZULENE = Path(__file__).parent / "files" / "azulene.xyz"


def test_inference_bonds_cov_radii():
    """
    Tries to infer connections from an .xyz using covalent radii.
    """
    chemgraph = cg.from_file(path_or_file=PATH_XYZ_AZULENE, fmt="xyz")

    assert chemgraph
    assert len(chemgraph.graph.edges) == 0

    chemgraph = chemgraph.infer_bonds(method="cov_radii")

    assert len(chemgraph.graph.edges) != 0


def test_inference_bonds_rdkit():
    """
    Tries to infer bonds from an .xyz using rdkit.
    """
    chemgraph = cg.from_file(path_or_file=PATH_XYZ_AZULENE, fmt="xyz")

    assert chemgraph
    assert len(chemgraph.graph.edges) == 0

    chemgraph = chemgraph.infer_bonds(method="rdkit")

    assert len(chemgraph.graph.edges) != 0
