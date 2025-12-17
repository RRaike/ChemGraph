from chemgraph.chemgraph import ChemGraph as cg
from chemgraph.fragment import detection

import rdkit.Chem

smiles = "c1ccccc1C=CC#C"
fragment_smarts = "C#C"
rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
rdkit_frag = rdkit.Chem.rdmolfiles.MolFromSmiles(fragment_smarts)

chemgraph = cg.from_file(rdkit_mol, fmt="mol")
cg_frag = cg.from_file(rdkit_frag, fmt="mol")

dict_fragments_rdkit = detection.detect_fragments_rdkit(
    chemgraph_or_graph=chemgraph, fragment_smarts=fragment_smarts
)

dict_fragments_nx = detection._detection_subgraph_(
    graph=chemgraph.graph, subgraph=cg_frag.graph
)


print(chemgraph)
