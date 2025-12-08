from .registry import register_geometry_parser
from ... import chemgraph
from ...utils import math, pathfinder
import networkx as nx


@register_geometry_parser("dihedrals")
def parse_dihedrals(chempgraph_or_graph: chemgraph.ChemGraph | nx.Graph):
    """
    Parses all dihedral angles in a ChemGraph or nx.Graph object.

    Args:
    -----
        chemgraph_or_graph: ChemGraph | nx.Graph

    Returns:
    --------
        list
    """
    g = chempgraph_or_graph
    if isinstance(g, chemgraph.ChemGraph):
        g = g.graph

    list_dihedral_angles = []

    list_paths = pathfinder._paths_finder_rev(
        g=g, n=3
    )  # All unique paths length 3 (= 4 nodes)

    for path in list_paths:
        ind_node_1 = path[0]
        ind_node_2 = path[1]
        ind_node_3 = path[2]
        ind_node_4 = path[3]

        node_1 = g.nodes[ind_node_1]
        node_2 = g.nodes[ind_node_2]
        node_3 = g.nodes[ind_node_3]
        node_4 = g.nodes[ind_node_4]

        dihedral_angle = math.dihedral_angle(
            pos_1=node_1["position"],
            pos_2=node_2["position"],
            pos_3=node_3["position"],
            pos_4=node_4["position"],
        )

        list_dihedral_angles.append(
            [(ind_node_1, ind_node_2, ind_node_3, ind_node_4), dihedral_angle]
        )

    return list_dihedral_angles
