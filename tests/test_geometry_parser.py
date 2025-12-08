from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path

PATH_XYZ_CYCLOHEXANE = Path(__file__).parent / "files" / "cyclohexane.xyz"


def test_geometry_parser():
    chemgraph = cg.from_file(path_or_file=PATH_XYZ_CYCLOHEXANE, fmt="xyz")
    chemgraph = chemgraph.infer_bonds(method="cov_radii")

    parsed_bonds = chemgraph.parse_geometry(geometry_parser="bonds")
    parsed_angles = chemgraph.parse_geometry(geometry_parser="angles")
    parsed_dihedrals = chemgraph.parse_geometry(geometry_parser="dihedrals")

    parsed_all = chemgraph.parse_geometry(
        geometry_parser=["bonds", "angles", "dihedrals"]
    )

    assert len(parsed_bonds["bonds"]) == len(chemgraph.graph.edges)
    assert len(parsed_angles["angles"]) > 1
    assert len(parsed_dihedrals["dihedrals"]) > 1
    assert len(parsed_all) == 3


#
# def test_bond_parser():
