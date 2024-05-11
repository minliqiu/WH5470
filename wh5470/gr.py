import numpy as np

G = 6.67430e-11  # Gravitational constant (N * m^2 / kg^2)
c = 299792458.0  # Speed of light (m/s)

def gr_modification(positions, velocities, masses):
    """
    Calculate the general relativistic modification to the acceleration.

    Args:
        positions (np.ndarray): Positions in Cartesian coordinates.
        velocities (np.ndarray): Velocities in Cartesian coordinates.
        masses (np.ndarray): Masses of bodies.

    Returns:
        np.ndarray: General relativistic modification to the acceleration.
    """
    num_bodies = len(masses)
    modification = np.zeros_like(positions)

    for i in range(num_bodies):
        for j in range(i + 1, num_bodies):
            r_ij = positions[j] - positions[i]
            r_ij_norm = np.linalg.norm(r_ij)
            v_ij = velocities[j] - velocities[i]
            v_ij_dot_r_ij = np.dot(v_ij, r_ij)

            modification[i] += G * masses[j] * (
                (4 * G * masses[j] / (c ** 2 * r_ij_norm)) * (
                    v_ij_dot_r_ij / r_ij_norm - 3 * (v_ij_dot_r_ij / r_ij_norm) ** 2 / (c ** 2)
                ) * r_ij / r_ij_norm
            )
            modification[j] -= G * masses[i] * (
                (4 * G * masses[i] / (c ** 2 * r_ij_norm)) * (
                    v_ij_dot_r_ij / r_ij_norm - 3 * (v_ij_dot_r_ij / r_ij_norm) ** 2 / (c ** 2)
                ) * r_ij / r_ij_norm
            )

    return modification