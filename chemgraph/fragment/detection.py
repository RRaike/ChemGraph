from .. import chemgraph
import networkx as nx
from typing import List

import rdkit.Chem


def detect_fragments_rdkit(
    chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph, fragment_smarts: str | List[str]
):
    """
    Detection of fragments in a ChemGraph or nx.Graph object using RDKit.

    Args:
    -----
        chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph
        fragment_smarts: str | List[str]
            SMARTS or list of SMARTS representations of the fragments to detect.

    Returns:
    --------
        dict(fragment_smarts_0: [(node_ind_0, ...),],
             fragment_smarts_1: ...,
             )
    """
    g = chemgraph_or_graph
    if isinstance(chemgraph_or_graph, nx.Graph):
        g = chemgraph.ChemGraph(name="g", graph=chemgraph_or_graph)

    assert isinstance(g, chemgraph.ChemGraph)

    mol = g.to_file(fmt="mol")
    dict_fragments = _detection_rdkit_(rdkit_mol=mol, fragment_smarts=fragment_smarts)

    return dict_fragments


def _detection_rdkit_(
    rdkit_mol: rdkit.Chem.rdchem.Mol, fragment_smarts: str | List[str]
):
    """
    Detection of fragments in a rdkit molecule.

    Args:
    -----
        rdkit_mol: rdkit.Chem.rdchem.Mol
            Representation of the molecule in question.
        fragment_smarts: str | List[str]
            SMARTS or list of SMARTS representations of the fragments to detect.

    Returns:
    --------
        dict(fragment_smarts_0: [(node_ind_0, ...),],
             fragment_smarts_1: ...,
             )
    """
    if isinstance(fragment_smarts, str):
        fragment_smarts = [
            fragment_smarts,
        ]

    dict_functional = dict()

    for frag in fragment_smarts:
        mol_fragment = rdkit.Chem.MolFromSmarts(frag)

        all_matches = rdkit_mol.GetSubstructMatches(
            mol_fragment
        )  # Subtructure search using RDkit.

        if len(all_matches) > 0:
            for match in all_matches:
                if frag not in dict_functional:
                    dict_functional[frag] = [
                        match,
                    ]
                else:
                    dict_functional[frag].append(match)

    return dict_functional


def _detection_subgraph_(graph: nx.Graph, subgraph: nx.Graph | List[nx.Graph]):
    """
    Detection of subgraphs in a nx.Graph.

    Args:
    -----
        graph: nx.Graph
            Graph wherein to look for pattern.
        subgraph: nx.Graph
            Subgraph to detect in the graph.

    Returns:
    --------

    """
    if not isinstance(subgraph, list):
        subgraph = [
            subgraph,
        ]

    # ------------------------------------------- #
    dict_subgraph = dict()

    for ind_subg, subg in enumerate(subgraph):
        gm = nx.isomorphism.GraphMatcher(
            graph,
            subg,
            node_match=nx.isomorphism.numerical_node_match("atom_number", 6),
            edge_match=nx.isomorphism.numerical_edge_match("bond_order", 1),
        )

        if gm.subgraph_is_isomorphic():
            if "name" in subg.graph:
                if subg.graph["name"] is not None:
                    key_name = subgraph.graph["name"]
                else:
                    key_name = ind_subg
            else:
                key_name = ind_subg

            matches = set()  # Set datatype prevents reverse graphs to be found. E.g. Node 0, 1 and Node 1, 0

            for subg_match in gm.subgraph_isomorphisms_iter():
                match = frozenset(set(subg_match.keys()))
                matches.add(match)

            if key_name not in dict_subgraph:
                dict_subgraph[key_name] = set()

            matches = [list(match) for match in matches]

            dict_subgraph[key_name] = matches

    return dict_subgraph
