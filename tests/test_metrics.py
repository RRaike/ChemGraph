from chemgraph.metrics import flexibility
from chemgraph.chemgraph import ChemGraph as cg
import rdkit.Chem


def test_flex_kier_alpha():
    """
    Tests Kier-Alpha flexibility metrics for graphs and chemgraphs.
    """
    smiles = "c1ccccc1C=CC#C"

    rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    chemgraph = cg.from_file(rdkit_mol, fmt="mol")

    for mode in ["a", "legacy"]:  # ['a', 'b', 'legacy]
        kier_alpha_cg = flexibility.kier_alpha(chemgraph_or_graph=chemgraph, mode=mode)
        kier_alpha_g = flexibility.kier_alpha(
            chemgraph_or_graph=chemgraph.graph, mode=mode
        )

        assert isinstance(kier_alpha_cg, float)
        assert isinstance(kier_alpha_g, float)
        assert kier_alpha_cg == kier_alpha_g


def test_shannon_entropy():
    """
    Tests Shannon entropy flexibility metrics for graphs and chemgraphs.
    """
    smiles = "c1ccccc1C=CC#C"

    rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    chemgraph = cg.from_file(rdkit_mol, fmt="mol")

    shannon_entropy_cg = flexibility.molecular_shannon_i(chemgraph_or_graph=chemgraph)
    shannon_entropy_g = flexibility.molecular_shannon_i(
        chemgraph_or_graph=chemgraph.graph
    )

    assert isinstance(shannon_entropy_cg, float)
    assert isinstance(shannon_entropy_g, float)
    assert shannon_entropy_cg == shannon_entropy_g


def test_kier_mkappa():
    """
    Tests Kier mkappa flexibility metrics for graphs and chemgraphs.
    """

    smiles = "c1ccccc1C=CC#C"

    rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    chemgraph = cg.from_file(rdkit_mol, fmt="mol")

    for m_opt in [0, 1, 2, 3]:
        for mode_opt in ["a", "legacy"]:
            for a_opt in [True, False]:
                kier_kappa_cg = flexibility.kier_mkappa(
                    chemgraph_or_graph=chemgraph, m=m_opt, alpha=a_opt, mode=mode_opt
                )

                kier_kappa_g = flexibility.kier_mkappa(
                    chemgraph_or_graph=chemgraph.graph,
                    m=m_opt,
                    alpha=a_opt,
                    mode=mode_opt,
                )

                assert isinstance(kier_kappa_cg, float)
                assert isinstance(kier_kappa_g, float)
                assert kier_kappa_cg == kier_kappa_g


def test_kier_phi():
    """
    Tests Kier Phi felexibility metrics for graphs and chemgraphs.
    """
    smiles = "c1ccccc1C=CC#C"

    rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    chemgraph = cg.from_file(rdkit_mol, fmt="mol")
    graph = chemgraph.graph

    kier_phi_cg = flexibility.kier_phi(
        chemgraph_or_graph=chemgraph, alpha=False, mode="a"
    )
    assert isinstance(kier_phi_cg, float)

    kier_phi_g = flexibility.kier_phi(chemgraph_or_graph=graph, alpha=False, mode="a")
    assert isinstance(kier_phi_g, float)

    assert kier_phi_cg == kier_phi_g
