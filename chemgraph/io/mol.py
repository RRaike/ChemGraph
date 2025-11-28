from .registry import register_reader, register_writer
import networkx as nx

import rdkit.Chem
from rdkit.Geometry import Point3D


RDKIT_TO_BO = {
    rdkit.Chem.rdchem.BondType.SINGLE: 1,
    rdkit.Chem.rdchem.BondType.DOUBLE: 2,
    rdkit.Chem.rdchem.BondType.TRIPLE: 3,
    rdkit.Chem.rdchem.BondType.AROMATIC: 1.5,
}

BO_TO_RDKIT = {v: k for k, v in RDKIT_TO_BO.items()}


@register_reader("mol")
def read_mol(mol: rdkit.Chem.rdchem.Mol, ind_conformer=0) -> dict:
    """
    Reads a mol rdkit.Chem.Mol object into a ChemGraph object.
    If the

    Args:
    -----
        mol: rdkit.Chem.rdchem.Mol()
            Molecule.


    """
    has_positions = False
    graph = nx.Graph()

    if (
        mol.GetNumConformers() != 0
    ):  # Checks if positional information is stored in the mol object.
        has_positions = True
        conf = list(mol.GetConformers())[ind_conformer]

    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    # atoms_kit = list(mol.GetAtoms())

    for atom in atoms:
        atom_number = atom.GetAtomicNum()
        position = None

        if has_positions:
            position = conf.GetAtomPosition(atom.GetIdx())

        graph.add_node(
            node_for_adding=atom.GetIdx(), atom_number=atom_number, position=position
        )

    for ind_bond, bond in enumerate(bonds):
        bond_order = RDKIT_TO_BO[bond.GetBondType()]

        graph.add_edge(
            u_of_edge=bond.GetBeginAtomIdx(),
            v_of_edge=bond.GetEndAtomIdx(),
            bond_order=bond_order,
        )

    return {"name": "from_mol", "graph": graph}


@register_writer("mol")
def write_mol(chemgraph) -> rdkit.Chem.rdchem.Mol:
    """
    Writes a ChemGraph object into a rdkit.Chem.Mol object.

    Args:
    -----
        chemgraph: ChemGraph
            ChemGraph object to be converted.

    Returns:
    --------
        mol: rdkit.Chem.rdchem.Mol
    """

    # ==== Create an editable molecule ==== #
    mol = rdkit.Chem.RWMol()
    has_positions = False

    for node, node_data in chemgraph.graph.nodes(data=True):
        atom = rdkit.Chem.Atom(node_data["atom_number"])
        mol.AddAtom(atom)

        if node_data["position"] is not None:
            has_positions = True

    for node_1, node_2, edge_data in chemgraph.graph.edges(data=True):
        bond_order = edge_data["bond_order"]
        rdkit_order = BO_TO_RDKIT[bond_order]
        mol.AddBond(node_1, node_2, rdkit_order)

    if has_positions:
        conf = rdkit.Chem.Conformer(len(chemgraph.graph.nodes()))

        for node, node_data in chemgraph.graph.nodes(data=True):
            position = node_data["position"]
            conf.SetAtomPosition(node, Point3D(position[0], position[1], position[2]))

        mol.AddConformer(conf)

    # Finalize molecule
    mol = mol.GetMol()
    rdkit.Chem.SanitizeMol(mol)

    return mol
