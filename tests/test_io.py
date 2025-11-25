import pytest

from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path

import rdkit.Chem

def test_non_existent_file():
    with pytest.raises(FileNotFoundError):
        path_xyz = Path(__file__).parent / "files" / "non_existent_file.xyz"
        cg.from_file(path_xyz)


def test_non_existent_reader():
    with pytest.raises(ValueError):
        path_xyz = path_xyz = Path(__file__).parent / "files" / "cyclohexane.xyz"
        cg.from_file(path_xyz, fmt="invalid_format")


def test_xyz_reader():
    path_xyz = Path(__file__).parent / "files" / "cyclohexane.xyz"
    chemgraph = cg.from_file(path_xyz)
    chemgraph_2 = cg.from_file(path_xyz, name="Test", fmt="xyz")

    assert chemgraph.name == path_xyz
    assert chemgraph_2.name == "Test"


def test_xyz_writer():
    path_xyz = Path(__file__).parent / "files" / "cyclohexane.xyz"
    path_write = path_xyz.parent / "cyclohexane_write.xyz"

    if path_write.exists():
        path_write.unlink()

    chemgraph = cg.from_file(path_xyz)
    chemgraph.to_file(path_write)
    chemgraph_2 = cg.from_file(path_write)

    assert chemgraph.graph.nodes(data=True) == chemgraph_2.graph.nodes(data=True)
    assert path_write.exists()
    path_write.unlink()



from rdkit.Chem import Draw
# from rdkit.Chem.Draw import IPythonConsole

mol = rdkit.Chem.rdmolfiles.MolFromSmiles('c1ccccc1C=CC#C')
# Draw.MolToFile(mol, 'test.png')
# print(mol)
chemgraph = cg.from_file(mol, fmt = 'mol')

