from .registry import register_geometry_parser
from ... import chemgraph
from ...utils import math, pathfinder
import networkx as nx


@register_geometry_parser("angles")
def parse_angles(chempgraph_or_graph: chemgraph.ChemGraph | nx.Graph):
    """
    Parses all bond angles in a ChemGraph or nx.Graph object.

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

    list_bond_angles = []
    list_paths = pathfinder._paths_finder_rev(
        g=g, n=2
    )  # All unique paths length 2 (= 3 nodes)

    for path in list_paths:
        ind_node_1 = path[0]
        ind_node_center = path[1]
        ind_node_2 = path[2]

        pos_1 = g.nodes[ind_node_1]["position"]
        pos_center = g.nodes[ind_node_center]["position"]
        pos_2 = g.nodes[ind_node_2]["position"]

        bond_angle = math.bond_angle(
            pos_center=pos_center,
            pos_1=pos_1,
            pos_2=pos_2,
        )

        list_bond_angles.append([(ind_node_1, ind_node_center, ind_node_2), bond_angle])

    return list_bond_angles
