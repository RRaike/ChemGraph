from .registry import register_reader, register_writer
import networkx as nx
import numpy as np

from pathlib import Path

from ..constants import periodic_table


@register_reader("xyz")
def read_xyz(path_xyz: str | Path) -> dict:
    """
    Reads a .xyz file into a name and a NetworkX graph.

    Args:
    -----
        path_xyz: str | Path
            Path to the .xyz file to read.

    Returns:
    --------
        dict
            {
            name: str | None.
                Name of the molecule. Takes the path as name.
            graph: nx.Graph
                Graph representation of the molecule.
                Does not infer bonds.
            }
    """
    with open(path_xyz, "r") as file:
        lines = file.readlines()

        comment = lines[1].strip()
        graph = nx.Graph()

        for ind_line, line in enumerate(lines[2:]):
            parts = line.split()

            if len(parts) > 4:
                raise ValueError("Invalid .xyz file format.")

            atom_type = parts[0]
            position = np.array(parts[1:]).astype(float)

            graph.add_node(
                node_for_adding=ind_line,
                atom_number=periodic_table.ATOMIC_NUM[atom_type],
                position=position,
            )

        graph.graph["description"] = comment

    return {"name": path_xyz, "graph": graph}


@register_writer("xyz")
def write_xyz(chemgraph, path: str | Path, precision="%22.15f"):
    """
    Writes a ChemGraph object to a .xyz file.

    Args:
    -----
        chemgraph: ChemGraph

        path: str | Path
            Path where to write the .xyz file.

    Returns:
    --------
        None
    """
    with open(path, "w+") as file:
        data_graph = chemgraph.graph.nodes(data=True)

        file.write(f"{len(data_graph)}\n")
        file.write(f"{chemgraph.graph.graph.get('description', '').strip()}\n")

        for node, data in data_graph:
            file.write(
                f"{periodic_table.ATOMIC_SYMBOLS[data['atom_number']]}"
                + f" {precision % data['position'][0]}"
                + f" {precision % data['position'][1]}"
                + f" {precision % data['position'][2]}\n"
            )

    return None
