import pytest
from chemgraph.chemgraph import ChemGraph as cg
from pathlib import Path

# def test_xyz_reader():
path_xyz = Path(__file__).parent / 'files' / 'cyclohexane.xyz'

a = cg.from_file(path_xyz)

b = 1