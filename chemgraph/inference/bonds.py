""" """

from .. import chemgraph
import networkx as nx
from rdkit.Chem import rdDetermineBonds
import ase.io
import ase.neighborlist

REGISTRY_INFERENCE_BONDS = dict()


def register_inference(name):
    def wrapper(func):
        REGISTRY_INFERENCE_BONDS[name] = func
        return func

    return wrapper


@register_inference("cov_radii")
def infer_bonds_cov_radii(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph,
):
    """
    Infers the bonds of a graph or ChemGraph using covalent radii powered by ASE.
    Removes all existing bonds before infering bonds.
    All bond orders are assumed to be 1.

    Args:
    -----
        chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph
            Representation of a molecule as a ChemGraph or a Graph.

    Returns:
    --------
        list

    """
    cg = chemgraph_or_graph
    if isinstance(chemgraph_or_graph, nx.Graph):
        cg = chemgraph.ChemGraph(name="g", graph=cg)

    assert isinstance(cg, chemgraph.ChemGraph)

    atoms = cg.to_file(fmt="atoms")
    cutoffs = ase.neighborlist.natural_cutoffs(atoms)
    nl = ase.neighborlist.NeighborList(
        cutoffs=cutoffs,
    )
    nl.update(atoms)

    edges = []

    for ind_atom_1 in range(len(atoms)):
        neigh, offset = nl.get_neighbors(ind_atom_1)

        for ind_atom_2 in neigh:
            if ind_atom_1 < ind_atom_2:
                edges.append((ind_atom_1, int(ind_atom_2), {"bond_order": 1}))

    return edges


@register_inference("rdkit")
def infer_bonds_rdkit(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph,
    charge: int = 0,
):
    """
    Infers the bonds of a graph using RDKIT.
    Removes all exisiting bonds before infering bonds.

    Args:
    -----
        chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph
            Representation of the molecule as ChemGraph or nx.Graph

        charge: (Optional) Int
            Default: 0
            Charge of the molecule.

    Returns:
    --------
        list
    """
    cg = chemgraph_or_graph
    if isinstance(cg, nx.Graph):
        cg = chemgraph.ChemGraph(name="graph", graph=cg)

    rdkit_mol = cg.to_file(fmt="mol")
    rdDetermineBonds.DetermineBonds(rdkit_mol, charge=charge)

    cg_rdkit = cg.from_file(path_or_file=rdkit_mol, fmt="mol")

    return list(cg_rdkit.graph.edges(data=True))
