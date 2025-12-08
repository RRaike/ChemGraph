import numpy as np

# ---------------------------------------------------------------------------------------------------------- #


def center_of_mass(positions: np.array, weights: np.array = None):
    """
    Calculates the center of mass of the positions with weights.
    If no weights are specified, they are assumed equivalent.

    Args:
    -----
        positions: np.array
            Positional arguments per row.

        weights: (Optional) np.array or list
            Default: None
            If specified, has to have one value per position.

    Returns:
    --------
        position_center_mass: np.array
    """
    if weights is None:  # If not specified, make row of 1s.
        weights = np.ones((len(positions), 1))

    elif isinstance(
        weights, list
    ):  # If specified, make sure to make it np.array for later computation.
        weights = np.array(weights)
        if len(positions) != len(weights):  # Check if the dimensions are good.
            raise ValueError("Positions and weights do not have same length!")

    positions_weighed = positions * weights

    position_center_mass = np.sum(positions_weighed, axis=0) / np.sum(weights)
    return position_center_mass


# ---------------------------------------------------------------------------------------------------------- #


def bond_length(pos_1: np.array, pos_2: np.array) -> float:
    """
    Returns the bond length between two atoms that are located at position_1 and position_2.

    Args:
    -----
        position_1: np.array
            Numpy array containing the position of atom 1.

        position_2: np.array
            Numpy array containing the position of atom 2.

    Returns:
    --------
        bond_length: Float
            Bond length between position_1 and position_2.
    """
    if pos_1.shape != pos_2.shape:
        raise TypeError("Positional vectors do not have the same shape.")

    bond_vector = np.subtract(pos_1, pos_2)  # Subtract the two positions.
    bond_length = np.linalg.norm(bond_vector)
    # bond_length = np.linalg.vector_norm(bond_vector)     # Take the vector norm of the difference.
    return bond_vector, bond_length


# ---------------------------------------------------------------------------------------------------------- #


def bond_angle(pos_center: np.array, pos_1: np.array, pos_2: np.array) -> float:
    """
    Returns the bond angle between the bonds of position_center, position_1 and position_center-position_2 in degrees.

    Args:
    -----
        position_center: np.array
            X, y, z position of the atom bound to the other atoms.

        position_1: np.array
            X, y, z position of an atom bound to the central atom.

        position_2: np.array
            X, y, z position of another atom bound to the central atom.

    Returns:
    --------
        bond_angle: float
    """
    if pos_center.shape != pos_1.shape or pos_center.shape != pos_2.shape:
        return TypeError("Vectors are not the same shape")
    elif np.array_equal(
        pos_1, pos_2
    ):  # Check if position 1 and 2 are the same. Edge case that breaks the calculation.
        return 0.0

    bond_vector_1, bond_length_1 = bond_length(pos_center, pos_1)
    bond_vector_2, bond_length_2 = bond_length(pos_center, pos_2)

    cosine_angle = np.dot(bond_vector_1, bond_vector_2) / (
        bond_length_1 * bond_length_2
    )
    bond_angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    bond_angle = np.rad2deg(bond_angle)

    return bond_angle


# ---------------------------------------------------------------------------------------------------------- #


def dihedral_angle(
    pos_1: np.array, pos_2: np.array, pos_3: np.array, pos_4: np.array
) -> float:
    """
    Returns the dihedral angle between the 4 atoms in degrees.

    Args:
    -----
        pos_1: np.array
            X, y, z coordinates of the first atom.

        pos_2: np.array
            X, y, z coordinates of the second atom.

        pos_3: np.array
            X, y, z coordinates of the third atom.

        pos_4: np.array
            X, y, z coordinates of the fourth atom.

    Returns:
    --------
        dihedral_angle: np.array
    """

    bond_1, bond_length_1 = bond_length(
        pos_1, pos_2
    )  # Sign is important here. Needs to be opposite of next two.
    bond_center, bond_length_center = bond_length(pos_3, pos_2)
    bond_2, bond_length_2 = bond_length(pos_3, pos_4)

    bond_center_unit = bond_center / bond_length_center  # Make unit vector.

    v = (
        bond_1 - np.dot(bond_1, bond_center_unit) * bond_center_unit
    )  # v = projection of bond_1 onto plane perpendicular to central bond
    #   = bond_1 minus component that aligns with central bond
    w = (
        bond_2 - np.dot(bond_2, bond_center_unit) * bond_center_unit
    )  # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1

    x = np.dot(v, w)  # v and w may not be normalized but that's fine since tan is y/x
    y = np.dot(
        np.cross(bond_center_unit, v), w
    )  # Angle between v and w in a plane is the dihedral angle
    dihedral_angle = np.arctan2(y, x)
    dihedral_angle = np.rad2deg(dihedral_angle)
    if dihedral_angle < 0:
        dihedral_angle += 360.0
    return dihedral_angle
