from chemgraph.chemgraph import ChemGraph as cg
from chemgraph.visualize import graph as vis_g
from pathlib import Path

PATH_XYZ_CYCLOHEXANE = Path(__file__).parent / "files" / "cyclohexane.xyz"
PATH_XYZ_AZULENE = Path(__file__).parent / "files" / "azulene.xyz"


chemgraph = cg.from_file(path_or_file=PATH_XYZ_AZULENE, fmt="xyz")

assert chemgraph
assert len(chemgraph.graph.edges) == 0

chemgraph = chemgraph.infer_bonds(method="rdkit")

fig = vis_g.draw_molecule_graph(graph=chemgraph.graph)

fig.show()


a = 1
