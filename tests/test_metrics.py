from chemgraph.metrics import flexibility
from chemgraph.chemgraph import ChemGraph as cg
import rdkit.Chem

smiles = "c1ccccc1C=CC#C"

rdkit_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
chemgraph = cg.from_file(rdkit_mol, fmt="mol")
a = flexibility.kier_alpha(chemgraph_or_graph=chemgraph, mode="legacy")

print(a)


shannon = flexibility.molecular_shannon_i(chemgraph_or_graph=chemgraph)
print(shannon)
