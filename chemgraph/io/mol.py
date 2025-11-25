from .registry import register_reader, register_writer
import networkx as nx
import numpy as np

import rdkit.Chem
from pathlib import Path

from .. import constants

RDKIT_TO_BOND_ORDERS = {
    rdkit.Chem.rdchem.BondType.SINGLE: 1,
    rdkit.Chem.rdchem.BondType.DOUBLE: 2,
    rdkit.Chem.rdchem.BondType.TRIPLE: 3,
    rdkit.Chem.rdchem.BondType.AROMATIC: 1.5
}

@register_reader("mol")
def read_mol(mol: rdkit.Chem.rdchem.Mol, ind_conformer = 0) -> dict:
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

    if mol.GetNumConformers() != 0: # Checks if positional information is stored in the mol object.
        has_positions = True
        conf = list(mol.GetConformers())[ind_conformer]
    
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    # atoms_kit = list(mol.GetAtoms())

    for ind_atom, atom in enumerate(atoms):
        atom_number = atom.GetAtomicNum()
        position = None

        if has_positions:  
            position = conf.GetAtomPosition(atom.GetIdx())

        graph.add_node(node_for_adding = atom.GetIdx(),
                       atom_number = atom_number,
                       position = position
                       )
        
    for ind_bond, bond in enumerate(bonds):
        bond_order = RDKIT_TO_BOND_ORDERS[bond.GetBondType()]

        graph.add_edge(
            u_of_edge  = bond.GetBeginAtomIdx(),
            v_of_edge  = bond.GetEndAtomIdx(),
            bond_order = bond_order 
        )
        
    a = 1

    return {
        'name': 'from_mol',
        'graph': graph
    }
