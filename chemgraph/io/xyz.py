from .registry import register_reader
import networkx as nx
import numpy as np

@register_reader('xyz')
def read_xyz(path_xyz):
    """
        Reads a .xyz file and returns a name and the atomic positions.
    """
    with open(path_xyz, 'r') as file:
        lines = file.readlines()

        comment = lines[1].strip()
        graph = nx.Graph()

        for ind_line, line in enumerate(lines[2:]):
            parts = line.split()

            if len(parts) > 4:
                raise ValueError("Invalid .xyz file format.")
            
            atom_type = parts[0]
            position  = np.array(parts[1:])

            graph.add_node(node_for_adding = ind_line, 
                           element = atom_type,
                           position = position
                           )

    return {
        'name': 'Test',
        'graph': graph
    }