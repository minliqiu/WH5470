import numpy as np

G = 6.67430e-11  # Gravitational constant (N m^2 / kg^2)

def cartesian_to_orbital(positions, velocities, masses):
    """
    Convert Cartesian coordinates to orbital elements.

    Args:
        positions (np.ndarray): Positions in Cartesian coordinates.
        velocities (np.ndarray): Velocities in Cartesian coordinates.
        masses (np.ndarray): Masses of bodies.

    Returns:
        tuple: Tuple containing orbital elements (a, e, i, Omega, omega, M) for each body.
    """
    orbital_elements = []
    for i in range(len(positions)):
        r = positions[i]
        v = velocities[i]
        total_mass = np.sum(masses)

        # Calculate angular momentum
        h = np.cross(r, v)
        h_norm = np.linalg.norm(h)

        # Calculate eccentricity vector
        k = np.cross(v, h) / G / total_mass - r / np.linalg.norm(r)
        e = np.linalg.norm(k)

        # Calculate semi-major axis
        a = 1 / (2 / np.linalg.norm(r) - np.linalg.norm(v) ** 2 / G / total_mass)

        # Calculate inclination
        i = np.arccos(h[2] / h_norm)

        # Calculate longitude of ascending node
        if h[0] >= 0 and h[1] >= 0:
            Omega = np.arctan2(h[1], h[0])
        elif h[0] < 0:
            Omega = np.arctan2(h[1], h[0]) + np.pi
        else:
            Omega = np.arctan2(h[1], h[0]) + 2 * np.pi

        # Calculate argument of periapsis
        if e > 1e-8:
            omega = np.arccos(np.dot(k, r) / e / np.linalg.norm(r))
            if np.cross(r, v)[2] < 0:
                omega = 2 * np.pi - omega
        else:
            omega = 0

        # Calculate mean anomaly
        E = np.arccos((a - np.linalg.norm(r)) / a / e)  # Eccentric anomaly
        if np.dot(r, v) < 0:
            E = 2 * np.pi - E
        M = E - e * np.sin(E)  # Mean anomaly

        orbital_elements.append((a, e, i, Omega, omega, M))

    return orbital_elements

def orbital_to_cartesian(orbital_elements, masses):
    """
    Convert orbital elements to Cartesian coordinates.

    Args:
        orbital_elements (tuple): Tuple containing orbital elements (a, e, i, Omega, omega, M) for each body.
        masses (np.ndarray): Masses of bodies.

    Returns:
        tuple: Tuple containing positions and velocities in Cartesian coordinates.
    """
    positions = []
    velocities = []
    total_mass = np.sum(masses)

    for a, e, i, Omega, omega, M in orbital_elements:
        # Calculate true anomaly
        E = M  # Initial guess for eccentric anomaly
        for _ in range(10):  # Newton-Raphson iteration
            E_new = E + (M - E + e * np.sin(E)) / (1 - e * np.cos(E))
            if np.abs(E_new - E) < 1e-8:
                break
            E = E_new
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

        # Calculate position and velocity
        r = a * (1 - e ** 2) / (1 + e * np.cos(nu))
        x = r * (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(i))
        y = r * (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(i))
        z = r * np.sin(omega + nu) * np.sin(i)
        position = np.array([x, y, z])

        v_r = np.sqrt(G * total_mass / a) * e * np.sin(nu) / (1 + e * np.cos(nu))
        v_theta = np.sqrt(G * total_mass / a) * (1 + e * np.cos(nu))
        v_x = (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(i)) * v_r \
              - (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(i)) * v_theta
        v_y = (np.sin(Omega) * np.cos(omega + nu) + np.cos(Omega) * np.sin(omega + nu) * np.cos(i)) * v_r \
              + (np.cos(Omega) * np.cos(omega + nu) - np.sin(Omega) * np.sin(omega + nu) * np.cos(i)) * v_theta
        v_z = np.sin(omega + nu) * np.sin(i) * v_r
        velocity = np.array([v_x, v_y, v_z])

        positions.append(position)
        velocities.append(velocity)

    return np.array(positions), np.array(velocities)