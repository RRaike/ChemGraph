from .registry import register_geometry_parser
from ... import chemgraph
from ...utils import math
import networkx as nx


@register_geometry_parser("bonds")
def parse_bonds(chempgraph_or_graph: chemgraph.ChemGraph | nx.Graph):
    """
    Parses all bonds in a ChemGraph or nx.Graph object.

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

    list_bond_length = []

    for ind_node_1, ind_node_2, edge_data in g.edges(data=True):
        pos_1 = g.nodes[ind_node_1]["position"]
        pos_2 = g.nodes[ind_node_2]["position"]

        _, bond_length = math.bond_length(pos_1=pos_1, pos_2=pos_2)

        list_bond_length.append([(ind_node_1, ind_node_2), bond_length])

    return list_bond_length
