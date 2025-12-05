import numpy as np
from itertools import compress
from collections import Counter
import warnings
import networkx as nx
import copy

from .. import chemgraph
from ..constants import periodic_table

# -------------------------------------------------------------------------------------- #

kier_radii = {
    (6, 4): 0.77,
    (6, 3): 0.67,
    (6, 2): 0.6,
    (7, 3): 0.74,
    (7, 2): 0.62,
    (7, 1): 0.55,
    (8, 2): 0.74,
    (8, 1): 0.62,
    (9, 1): 0.72,
    (15, 5): 1.1,
    (15, 4): 1.1,
    (15, 3): 1.1,
    (15, 2): 1.0,
    (15, 1): 0.95,  # Taken from P2  ISBN 978-0442233945
    (16, 4): 1.04,
    (16, 3): 1.04,
    (16, 2): 1.04,
    (16, 1): 0.94,
    (17, 1): 0.99,
    (35, 1): 1.14,
    (53, 1): 1.33,
}

# -------------------------------------------------------------------------------------- #


def kier_alpha(
    chemgraph_or_graph: nx.Graph | chemgraph.ChemGraph,
    radii: dict = periodic_table.COVALENT_RADII,
    mode: str = "a",
) -> float:
    """
    Compute the alpha correction factor for Kier indices based on atomic radii.

    Args:
    -----
        chemgraph_or_graph: ChemGraph | networkx.Graph
            Molecular graph.
        radii: dict
            Dictionary of covalent radii.
        mode: str
            Mode of alpha computation ('a', 'b', 'legacy').

    Returns:
        float: Alpha correction value.
    """
    g = copy.deepcopy(chemgraph_or_graph)
    if isinstance(chemgraph_or_graph, chemgraph.ChemGraph):
        g = chemgraph_or_graph.graph

    has_H = any(data.get("atom_number") == 1 for _, data in g.nodes(data=True))

    if has_H:
        indx_H = [
            ind_n
            for ind_n, _ in g.nodes(data=True)
            if g.nodes[ind_n["atom_number"]] == 1
        ]
        g.remove_nodes_from(indx_H)

    k_alpha = 0

    if mode == "a":
        for ind_note, node_data in g.nodes(data=True):
            # if atom_symbol == 'X':
            #     atom_symbol = 'C'
            atom_number = node_data["atom_number"]
            k_alpha += (radii[atom_number] / radii[6]) - 1
    elif mode == "b":
        for ind_node_1, ind_node_2, edge in g.edges(data=True):
            position_1 = g.nodes[ind_node_1]["position"]
            position_2 = g.nodes[ind_node_2]["position"]
            bond_length = np.linalg.norm(np.array(position_1) - np.array(position_2))
            k_alpha += (bond_length / 1.535) - 1  # Taken from Wiki, Sp3 C-C bond
    elif mode == "legacy":
        warnings.warn("Legacy mode. Use only for uncharged/non-radical molecules.")
        for ind_node, node_data in g.nodes(data=True):
            atom_number = node_data["atom_number"]
            # if atom_symbol == 'X':
            #     atom_symbol = 'C'

            if atom_number != 1:
                hybridization = g.degree[ind_node]

                if (atom_number, hybridization) in kier_radii:
                    kier_radius = kier_radii[(atom_number, hybridization)]

                if (atom_number, hybridization) not in kier_radii:
                    warnings.warn(
                        f"Atomic number '{atom_number}' not tabulated. Using sp3 carbon."
                    )
                    kier_radius = kier_radii[(6, 4)]

                k_alpha += (kier_radius / kier_radii[(6, 4)]) - 1

    else:
        raise NotImplementedError(f"No mode '{mode}'.")

    return k_alpha


# -------------------------------------------------------------------------------------- #


def molecular_shannon_i(chemgraph_or_graph: nx.Graph | chemgraph.ChemGraph) -> float:
    """
    Compute the Shannon entropy of the atom types in a graph.
    Calculates the Shannon entropy of the molecule without taking Hydrogens into account.

    Args:
    -----
        chemgraph_or_graph: ChemGraph | nx.Graph
            Molecular graph.

    Returns:
        float: Shannon entropy value.
    """
    # Unpack if ChemGraph
    g = copy.deepcopy(chemgraph_or_graph)
    if isinstance(chemgraph_or_graph, chemgraph.ChemGraph):
        g = chemgraph_or_graph.graph

    has_H = any(data.get("atom_number") == 1 for _, data in g.nodes(data=True))

    if has_H:
        indx_H = [
            ind_n
            for ind_n, _ in g.nodes(data=True)
            if g.nodes[ind_n["atom_number"]] == 1
        ]
        g.remove_nodes_from(indx_H)

    num_atoms = len(g.nodes())
    atom_types = []

    for ind_node, data_node in g.nodes(data=True):
        # Calculate connectivity of atom.
        # atom_connectivity = atom.atom_connectivity(graph = g, ind_atom = na)
        atom_connectivity = g.degree[ind_node]
        atom_number = data_node["atom_number"]

        # --------------------------- #
        # Get atom number.
        # atom_symbol = g.nodes[na]['atom_symbol']

        # if atom_symbol == 'X':
        #     atom_symbol = 'C'
        # atom_number = atomic_numbers[atom_symbol]

        # --------------------------- #

        atom_types.append(tuple([atom_number, atom_connectivity]))

    freqs = dict([(k, v / num_atoms) for k, v in Counter(atom_types).items()])

    i = 0
    for _, rho_i in freqs.items():
        i -= rho_i * np.log10(rho_i)

    return i


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


def kier_mkappa(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph,
    m: int,
    alpha: bool = False,
    mode: str = "a",
):
    """
    Compute the m-th order Kier kappa shape index.

    Args:
    -----
        chemgraph_or_graph: ChemGraph | nx.Graph
            Molecular graph.
        m: int
            Order of kappa index (0, 1, 2, or 3).
        alpha: bool
            Whether to include alpha correction.
        mode: str
            Mode for alpha correction ('a', 'b', 'legacy').

    Returns:
    --------
        float: m-th order kappa shape index.
    """
    # Unpack if ChemGraph
    g = copy.deepcopy(chemgraph_or_graph)
    if isinstance(chemgraph_or_graph, chemgraph.ChemGraph):
        g = chemgraph_or_graph.graph

    if mode == "legacy":
        indx_H = [
            ind_n
            for ind_n, _ in g.nodes(data=True)
            if g.nodes[ind_n]["atom_number"] == 1
        ]
        g.remove_nodes_from(indx_H)

    if alpha:
        alf = kier_alpha(g, mode=mode)
    else:
        alf = 0

    # if mode == "legacy":
    #     g = no_hydrogen(g)
    # a = g.number_of_nodes()
    # num_atoms, num_H = atom.number_of_atoms(graph = g)          # Number of atoms in the graph. If H is implicit, it will not count those in that number.
    num_atoms = len(g.nodes())

    if m == 0:
        return molecular_shannon_i(g) * num_atoms

    elif m == 1:
        num = (num_atoms + alf - 0) * (num_atoms + alf - 1) ** 2
        p = g.number_of_edges()

    elif m == 2:
        num = (num_atoms + alf - 1) * (num_atoms + alf - 2) ** 2
        unique_paths = _paths_finder_rev(g, m)
        p = len(unique_paths)

    elif m == 3:
        assert num_atoms > 2, f"Needs at least 3 atoms, got '{num_atoms}'."
        # TODO
        warnings.warn("3K may not works with cyclopropanes.")
        if num_atoms % 2 == 0:
            num = (num_atoms + alf - 3) * ((num_atoms + alf - 2) ** 2)
        else:
            num = (num_atoms + alf - 1) * ((num_atoms + alf - 3) ** 2)
        unique_paths = _paths_finder_rev(g, m)
        p = len(unique_paths)

    else:
        raise NotImplementedError("Invalid 'm', '{m}'.")

    return num / ((p + alf) ** 2)


# -------------------------------------------------------------------------------------- #


def kier_phi(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph,
    alpha: bool = False,
    mode: str = "a",
) -> float:
    """
    Compute the Kier phi descriptor of a molecule.
    Quantify molecular shape and branching in a size-normalized way.
    Higher kier_phi means a molecule is likely more branched and/or complex.
    Lower kier_phi means a molecule is more linear or simple.

    Args:
    -----
        chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph
            Molecular graph.
        alpha: bool
            Whether to include alpha correction.
        mode: str
            Mode for alpha correction ('a', 'b', 'legacy').

    Returns:
    --------
        float: Kier phi descriptor.
    """
    # Unpack if ChemGraph
    g = copy.deepcopy(chemgraph_or_graph)
    if isinstance(chemgraph_or_graph, chemgraph.ChemGraph):
        g = chemgraph_or_graph.graph

    if mode == "legacy":
        indx_H = [
            ind_n
            for ind_n, _ in g.nodes(data=True)
            if g.nodes[ind_n["atom_number"]] == 1
        ]
        g.remove_nodes_from(indx_H)

    num_of_atoms = len(g.nodes())

    # # if mode == "legacy":
    # #     num_a = no_hydrogen(g).number_of_nodes()
    # # else:
    # #     num_a = g.number_of_nodes()
    # # num_a = g.number_of_nodes()
    # num_of_atoms = atom.number_of_atoms(graph = g)
    return (
        kier_mkappa(g, 1, alpha=alpha, mode=mode)
        * kier_mkappa(g, 2, alpha=alpha, mode=mode)
        / num_of_atoms
    )


# -------------------------------------------------------------------------------------- #


def crest_flex(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph, bo_label="bond_order"
) -> float:
    """
    Compute the flexibility score of a molecule from its graph.

    The flexibility is computed as:
        Flexibility = sqrt((1/m) * sum(val^2 for all bonds))

    where 'val' for each bond depends on branching, ring membership, bond order, and hybridization.

    Args:
    -----
        chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph
            Molecular graph.
        bo_label: str
            Edge attribute name representing bond order. Default is 'bond_order'.

    Returns:
    --------
        float
    """
    g = copy.deepcopy(chemgraph_or_graph)
    if isinstance(chemgraph_or_graph, chemgraph.ChemGraph):
        g = chemgraph_or_graph.graph

    has_H = any(data.get("atom_number") == 1 for _, data in g.nodes(data=True))

    if has_H:
        indx_H = [
            ind_n
            for ind_n, _ in g.nodes(data=True)
            if g.nodes[ind_n["atom_number"]] == 1
        ]
        g.remove_nodes_from(indx_H)

    # if g.graph["graph_type"] != "kenogram":
    #     warnings.warn("Per definition 'crest_flex' works on kenograms.")
    if not nx.get_edge_attributes(g, bo_label) or bo_label == "":
        nx.set_edge_attributes(g, 0.0, bo_label)
        warnings.warn("'crest_flex' No bond orders found.")

    cycles = nx.cycle_basis(g)
    # cycles = structures.all_rings(graph = g, length_haptic_cycle = 3)     # Accounts for ghost nodes and haptic bonds.
    edges = list(g.edges())

    m = len(edges)
    av2 = 0.0
    for ind_node_1, ind_node_2, data_edge in g.edges(data=True):
        # bond_order = g.edges[edge][bo_label]
        bond_order = data_edge[bo_label]

        if bond_order != 0:
            # assert len(edge) == 2

            # ------------------------------------------------ #
            # Haptic bonds add rigidity and act in that fashion
            # as if they are single bonds. Accounting for this.

            # if g.edges[edge][bo_label] == 'h':
            #     bond_order = 1

            # ------------------------------------------------ #
            # cns = np.array([atom.atom_connectivity(graph = g, ind_atom = ind_atom) for ])
            # cns = np.array([len(list(g.neighbors(n))) for n in edge])

            cns = [g.degree[ind_node_1], g.degree[ind_node_2]]
            # for ind_atom in edge:
            #     connectivity = atom.atom_connectivity(graph = g, ind_atom = ind_atom)
            #     cns.append(connectivity)

            cns = np.array(cns)
            # ------------------------------------------------ #

            hybf = 1.0

            hybf *= (
                0.5 if g.nodes[ind_node_1]["atom_number"] == 6 and cns[0] < 4 else 1.0
            )
            hybf *= (
                0.5 if g.nodes[ind_node_2]["atom_number"] == 6 and cns[1] < 4 else 1.0
            )

            # for ind_edge in range(len(edge)):
            #     node = g.nodes[edge[ind_edge]]
            #     atom_number = atomic_numbers[node['atom_symbol']]

            #     if atom_number == 6 and cns[ind_edge] < 4:
            #         hybf*= 0.5

            # ------------------------------------------------ #

            doublef = 1.0 - np.exp(-4.0 * (bond_order - 2.0) ** 6)
            branch = 2.0 / np.sqrt(np.prod(cns))

            iring = [np.isin([ind_node_1, ind_node_2], x).any() for x in cycles]

            try:
                k = min([len(x) for x in list(compress(cycles, iring))])
            except ValueError:
                k = 0

            ringf = 0.5 * (1.0 - np.exp(-0.06 * k)) if k > 0 else 1.0
            val = branch * ringf * doublef * hybf
            av2 += val**2

    av2 = np.sqrt(av2 / m) if m > 0 else av2
    return av2


# -------------------------------------------------------------------------------------- #
