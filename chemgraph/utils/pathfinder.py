import networkx as nx

# -------------------------------------------------------------------------------------- #


def recu_path(
    g: nx.Graph, na: int, n: int, sub_paths: list = [], path: list | None = None
):
    """
    Recursive helper function to find all unique simple paths of length n starting from node na.

    Args:
    -----
        g: nx.Graph
        na: Int
            Node index
        n: Int
            Length of path remaining
        sub_paths: List of Lists
            Default: []
            Initiation for recursively finding sub paths.
        path: List | None
            Default: None
            Initiation for initial step.

    Returns:
    --------
        sub_paths: List of list
    """
    if path is None:
        path = [na]
    for neighbor in g.neighbors(path[-1]):
        if neighbor not in path:
            edge_attributes = g.get_edge_data(path[-1], neighbor)
            if n - 1 == 0 and edge_attributes["bond_order"] != 0:
                if neighbor > na:
                    sub_paths += [path + [neighbor]]
            else:
                recu_path(g, na, n - 1, sub_paths, path + [neighbor])
    return sub_paths


def _paths_finder_rev(g: nx.Graph, n: int):
    """
    Find all unique simple paths of length n in the graph.

    Args:
    -----
        g: networkx.Graph
            Molecular graph.
        n: int
            Path length.

    Returns:
    --------
        list: List of paths (each path is a list of node indices).
    """
    paths = []
    if n < 1:
        for na in g.nodes():
            paths.append([na])
    else:
        for na in g.nodes():
            paths.extend(recu_path(g, na, n, []))
    return paths


# -------------------------------------------------------------------------------------- #
