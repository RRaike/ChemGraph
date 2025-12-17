""" """

import networkx as nx
from ..constants import colors, periodic_table
import matplotlib.pyplot as plt


def draw_molecule_graph(
    graph,
    show_index=True,
):
    """
    Creates a matplotlib figure object of the molecular graph.

    Args:
    -----
        graph: nx.Graph()
            Graph representation of the molecule.
            Must contain 'symbol' keyword in its nodes and 'bond_order' in its edges.
        show_index: (Optional) Boolean
            Default: True
            If True, atom symbols will in subscript show their index in the graph.

    Returns:
    ---------
        fig: matplotlib.figure.Figure
            Figure object of the molecular graph.
    """
    # colors = conventions.CPK_COLOR
    styles_edges = {1: "solid", 2: "dashed", 3: "dotted", 1.5: (0, (1, 1))}
    list_colors = []
    dict_nodes = {}
    list_edges = []

    positions = nx.kamada_kawai_layout(graph, center=[0, 1])

    for ind_node, data_node in graph.nodes(data=True):
        # node_attributes = graph.nodes[node]
        # symbol = node_attributes['atom_symbol']
        atom_num = data_node["atom_number"]
        # symbol = symbol.replace('*', 'X')

        color = colors.CPK_COLOR_NUMBERS[atom_num]
        color_normalized = tuple([c / 256 for c in color])

        list_colors.append(color_normalized)
        dict_nodes[ind_node] = (
            periodic_table.ATOMIC_SYMBOLS[atom_num] + "$_{" + str(ind_node) + "}$"
        )

    for _, _, edge_data in graph.edges(data=True):
        bond_order = edge_data["bond_order"]

        style = "solid"

        if bond_order in styles_edges:
            style = styles_edges[bond_order]

        list_edges.append(style)
        # list_edges.append(edge_atrributes['bond_order'])

    fig = plt.figure()

    nx.draw(
        graph,
        node_size=350,
        labels=dict_nodes,
        edgecolors="black",
        node_color=list_colors,
        style=list_edges,
        pos=positions,
        alpha=0.75,
    )

    return fig
