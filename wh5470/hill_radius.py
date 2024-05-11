import numpy as np

G = 6.67430e-11  # Gravitational constant (N * m^2 / kg^2)

def calculate_hill_radius(m1, m2, a):
    """
    Calculate the Hill radius for a body orbiting another body.

    Args:
        m1 (float): Mass of the primary body.
        m2 (float): Mass of the secondary body.
        a (float): Semi-major axis of the orbit.

    Returns:
        float: Hill radius.
    """
    return a * ( m2/(3*(m1+m2)) )**(1 / 3)