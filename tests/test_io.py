import pytest

from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path

import rdkit.Chem

import ase.io


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

    assert path_write.exists()
    path_write.unlink()


def test_io_mol():
    smiles = "c1ccccc1C=CC#C"
    rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    smiles_restructured = rdkit.Chem.MolToSmiles(
        rdkit_mol
    )  # RDkit can reshuffle a SMILES

    chemgraph = cg.from_file(rdkit_mol, fmt="mol")

    assert len(chemgraph.graph.nodes()) == rdkit_mol.GetNumAtoms()

    mol_obj = chemgraph.to_file(fmt="mol")
    smiles_cg = rdkit.Chem.MolToSmiles(mol_obj)

    assert smiles_restructured == smiles_cg


def test_io_atoms():
    PATH_XYZ_AZULENE = Path(__file__).parent / "files" / "azulene.xyz"
    atoms = ase.io.read(PATH_XYZ_AZULENE, index=-1)

    chemgraph = cg.from_file(atoms, fmt="atoms")
    atoms_cg = chemgraph.to_file(fmt="atoms")

    assert atoms == atoms_cg
