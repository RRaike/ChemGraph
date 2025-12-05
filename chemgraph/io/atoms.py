"""
IO from and to ase.Atoms objects.
"""

from .registry import register_reader, register_writer
from .. import chemgraph
import ase
import ase.io
import networkx as nx


@register_reader("atoms")
def read_atoms(ase_atoms: ase.Atoms) -> dict:
    """
    Reads an ASE.Atoms object into a name and a Networkx.Graph

    Args:
    -----
        ase_atoms: ase.Atoms
            Representation of Molecule

    Returns:
    --------
        dict: {'name': str, 'graph: nx.Graph}
    """
    name = ase_atoms.get_chemical_formula()
    graph = nx.Graph()

    atom_numbers = ase_atoms.get_atomic_numbers()
    positions = ase_atoms.get_positions()

    for ind_number, atom_number in enumerate(atom_numbers):
        graph.add_node(
            node_for_adding=ind_number,
            atom_number=int(atom_number),
            position=positions[ind_number, :],
        )

    return {"name": name, "graph": graph}


@register_writer("atoms")
def write_atoms(chemgraph: chemgraph.ChemGraph) -> ase.Atoms:
    """
    Writes a ChemGraph into a ASE.Atoms object.

    Args:
    -----
        chemgraph: chemgraph.ChemGraph
            ChemGraph representation of a molecule.

    Returns:
    --------
        ase.Atoms
    """
    atom_numbers = []
    positions = []

    for _, data in chemgraph.graph.nodes(data=True):
        atom_numbers.append(data["atom_number"])
        positions.append(data["position"])

    atoms = ase.Atoms(numbers=atom_numbers, positions=positions)
    return atoms
