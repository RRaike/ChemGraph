"""
File adds elements to an existing gmol object
"""


# def infer_geometries_c2m(
#     chemgraph_or_graph: chemgraph.ChemGraph | nx.Graph,
# ):
#     """
#     Iterates over all atoms in the gmol object and assigns a molecular geometry to every atom in the object.
#     Geometries are stored in the atom object.
#     Alters input gmol object

#     Args:
#     -----
#         gmol: gmol object

#     Returns:
#     --------
#         gmol: gmol object
#     """
#     pos_all = (
#         gmol.coord
#     )  # Faster to cycle through a list than call attributes from objects.
#     labels_all = gmol.labels  # Idem

#     for ind_atom, atom in enumerate(gmol.atoms):
#         list_positions = []
#         list_labels = []

#         list_labels.append(atom.label)
#         list_positions.append(atom.coord)

#         indx_neighbors = atom.adjacency

#         for ind_neighbor in indx_neighbors:
#             label_neighbor = labels_all[ind_neighbor]
#             pos_neighbor = pos_all[ind_neighbor]

#             list_labels.append(label_neighbor)
#             list_positions.append(pos_neighbor)

#         posgeom_dev = coordination_sphere.shape_measure(list_labels, list_positions)
#         coordination_geometry = min(posgeom_dev, key=posgeom_dev.get)

#         atom.geometry = coordination_geometry
#         atom.geometry_dev = posgeom_dev[coordination_geometry]

#     return gmol
