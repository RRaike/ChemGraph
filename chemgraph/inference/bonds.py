""" """

from .. import chemgraph
import networkx as nx
from rdkit.Chem import rdDetermineBonds


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

    Args:
    -----
    """

    return


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
    """
    cg = chemgraph_or_graph
    if isinstance(cg, nx.Graph):
        cg = chemgraph.ChemGraph(name="graph", graph=cg)

    rdkit_mol = cg.to_file(fmt="mol")
    rdDetermineBonds.DetermineBonds(rdkit_mol, charge=charge)

    cg = chemgraph_or_graph.from_file(path_or_file=rdkit_mol, fmt="mol")

    return list(cg.graph.edges(data=True))
