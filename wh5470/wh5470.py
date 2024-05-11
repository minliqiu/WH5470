import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from orbital_elements import cartesian_to_orbital, orbital_to_cartesian
from hill_radius import calculate_hill_radius
from gr import gr_modification

G = 6.67430e-11  # Gravitational constant (N m^2 / kg^2)

class wh5470:
    def __init__(self, masses, orbital_elements=None, positions=None, velocities=None, time_step=0.1, total_time=100., gr=False):
        """
        Initialize the N-body simulator.

        Args:
            masses (np.ndarray): Array of masses for each body.
            orbital_elements (tuple, optional): Tuple containing orbital elements (a, e, i, Omega, omega, M) for each body.
            positions (np.ndarray, optional): Initial positions of bodies in Cartesian coordinates.
            velocities (np.ndarray, optional): Initial velocities of bodies in Cartesian coordinates.
            time_step (float): Time step for the integration.
            total_time (float): Total time to simulate.
        """
        self.masses = masses
        self.positions = positions
        self.velocities = velocities
        self.time_step = time_step
        self.total_time = total_time
        self.num_bodies = len(masses)
        self.gr = gr
        
        if orbital_elements is not None:
            self.positions, self.velocities = orbital_to_cartesian(orbital_elements, masses)
        else:
            self.positions = positions
            self.velocities = velocities

    def cartesian_to_jacobi(self, positions, velocities):
        """
        Convert Cartesian coordinates to Jacobi coordinates.

        Args:
            positions (np.ndarray): Positions in Cartesian coordinates.
            velocities (np.ndarray): Velocities in Cartesian coordinates.

        Returns:
            tuple: Tuple containing Jacobi positions and velocities.
        """
        jacobi_positions = np.zeros_like(positions)
        jacobi_velocities = np.zeros_like(velocities)

        com_position = np.sum(self.masses[:, None] * positions, axis=0) / np.sum(self.masses)
        com_velocity = np.sum(self.masses[:, None] * velocities, axis=0) / np.sum(self.masses)

        for i in range(self.num_bodies):
            jacobi_positions[i] = positions[i] - com_position
            jacobi_velocities[i] = velocities[i] - com_velocity

        return jacobi_positions, jacobi_velocities

    def jacobi_to_cartesian(self, jacobi_positions, jacobi_velocities):
        """
        Convert Jacobi coordinates to Cartesian coordinates.

        Args:
            jacobi_positions (np.ndarray): Positions in Jacobi coordinates.
            jacobi_velocities (np.ndarray): Velocities in Jacobi coordinates.

        Returns:
            tuple: Tuple containing Cartesian positions and velocities.
        """
        cartesian_positions = np.zeros_like(jacobi_positions)
        cartesian_velocities = np.zeros_like(jacobi_velocities)

        com_position = np.sum(self.masses[:, None] * jacobi_positions, axis=0) / np.sum(self.masses)
        com_velocity = np.sum(self.masses[:, None] * jacobi_velocities, axis=0) / np.sum(self.masses)

        for i in range(self.num_bodies):
            cartesian_positions[i] = jacobi_positions[i] + com_position
            cartesian_velocities[i] = jacobi_velocities[i] + com_velocity

        return cartesian_positions, cartesian_velocities

    def acceleration(self, positions):
        """
        Calculate the acceleration for each body.

        Args:
            positions (np.ndarray): Positions in Jacobi coordinates.

        Returns:
            np.ndarray: Accelerations for each body.
        """
        accelerations = np.zeros_like(positions)
        for i in range(self.num_bodies):
            for j in range(self.num_bodies):
                if i != j:
                    r_ij = positions[j] - positions[i]  # vector from body i to body j
                    dist_ij = np.linalg.norm(r_ij)  # distance between bodies i and j
                    accelerations[i] += G * self.masses[j] * r_ij / dist_ij**3
        
        
#         for i in range(self.num_bodies):
#             for j in range(i + 1, self.num_bodies):
#                 r_ij = positions[j] - positions[i]
#                 r_ij_norm = np.linalg.norm(r_ij)
#                 accelerations[i] -= G * self.masses[j] * r_ij / (r_ij_norm ** 3)
#                 accelerations[j] += G * self.masses[i] * r_ij / (r_ij_norm ** 3)
        return accelerations
    
    def step(self):
        '''
        split the step in Cartesian coordinates
        '''
        self.velocities += 0.5 * self.acceleration(self.positions) * self.time_step
        self.positions += self.velocities * self.time_step
        self.velocities += 0.5 * self.acceleration(self.positions) * self.time_step
        position_cartisian, velocity_cartisian = self.positions, self.velocities
        return position_cartisian, velocity_cartisian

    def integrate(self):
        """
        Perform the Wisdom-Holman integrator (kick-drift-kick algorithm).
        """
        time = 0.0
        jacobi_positions, jacobi_velocities = self.cartesian_to_jacobi(self.positions, self.velocities)
        trajectory = [self.positions.copy()]
        velocity = [self.velocities.copy()]
        orbital_elements_array = [cartesian_to_orbital(self.positions, self.velocities, self.masses).copy()]
        energy, angular_momentum = [], []
        trajectory_jacobi, velocity_jacobi = [], []

        while time < self.total_time:
            orbital_elements = cartesian_to_orbital(self.positions, self.velocities, self.masses)
            # Check for Hill radius violation
            for i in range(1, self.num_bodies):
                for j in range(i + 1, self.num_bodies):
                    r_ij = np.linalg.norm(self.positions[j] - self.positions[i])
                    a = orbital_elements[i][0] # Semi-major axis for body i
                    hill_radius = calculate_hill_radius(self.masses[i], self.masses[j], a)
                    if r_ij < hill_radius:
                        print(f"Body {j + 1} entered the Hill radius of Body {i + 1} at time {time:.2f}")
                        return
            position_cart, velocity_cart = self.step()
            # First kick
            jacobi_accelerations = self.acceleration(jacobi_positions)
            if self.gr:
                cartesian_accelerations = self.jacobi_to_cartesian(jacobi_positions, jacobi_velocities)[1]
                gr_modification_cartesian = gr_modification(self.positions, self.velocities, self.masses)
                gr_modification_jacobi = self.cartesian_to_jacobi(self.positions, gr_modification_cartesian)[1]
                jacobi_accelerations += gr_modification_jacobi
            jacobi_velocities += 0.5 * self.time_step * jacobi_accelerations

            # Drift
            jacobi_positions += self.time_step * jacobi_velocities

            # Second kick
            jacobi_accelerations = self.acceleration(jacobi_positions)
            if self.gr:
                cartesian_accelerations = self.jacobi_to_cartesian(jacobi_positions, jacobi_velocities)[1]
                gr_modification_cartesian = gr_modification(self.positions, self.velocities, self.masses)
                gr_modification_jacobi = self.cartesian_to_jacobi(self.positions, gr_modification_cartesian)[1]
                jacobi_accelerations += gr_modification_jacobi
            jacobi_velocities += 0.5 * self.time_step * jacobi_accelerations

            time += self.time_step

            self.positions, self.velocities = position_cart, velocity_cart # self.jacobi_to_cartesian(jacobi_positions, jacobi_velocities)
            trajectory.append(self.positions.copy())
            velocity.append(self.velocities.copy())
            energy.append(self.calculate_total_energy())
            angular_momentum.append(self.calculate_total_angular_momentum())
            orbital_elements_array.append(cartesian_to_orbital(self.positions, self.velocities, self.masses))
            p_jacobi, v_jacobi = self.cartesian_to_jacobi(self.positions, self.velocities)
            trajectory_jacobi.append(p_jacobi)
            velocity_jacobi.append(v_jacobi)
        return trajectory, velocity, energy, angular_momentum, orbital_elements_array, trajectory_jacobi, velocity_jacobi


    def plot_trajectory(self, trajectory):
            """
            Plot the trajectory of the N-body simulation.

            Args:
                trajectory (np.ndarray): Array containing the positions of bodies at each time step.
            """
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            
            iter_array = range(len(trajectory))

            for i in range(self.num_bodies):
                x = [trajectory[j][i][0] for j in iter_array]
                y = [trajectory[j][i][1] for j in iter_array]
                z = [trajectory[j][i][2] for j in iter_array]
                ax.plot(x, y, z, label=f'Body {i + 1}')

            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title('N-Body Simulation Trajectory')
            ax.legend()

            plt.show()

    def calculate_total_energy(self):
        total_kinetic_energy = 0.5 * np.sum(self.masses * np.sum(self.velocities**2, axis=1))
        total_potential_energy = 0
        for i in range(self.num_bodies):
            for j in range(i + 1, self.num_bodies):
                distance = np.linalg.norm(self.positions[j] - self.positions[i])
                total_potential_energy -= G * self.masses[i] * self.masses[j] / distance
        return total_kinetic_energy + total_potential_energy

    def calculate_total_angular_momentum(self):
        total_angular_momentum = np.zeros(3)
        for i in range(self.num_bodies):
            total_angular_momentum += np.cross(self.positions[i], self.masses[i] * self.velocities[i])
        return total_angular_momentum

    def plot_energy_and_momentum(self, energy, angular_momentum):
        fig, axs = plt.subplots(2, 1, figsize=(8, 6))
        axs[0].set_ylim([1.1*energy[0], 0.9*energy[0]])
        axs[0].plot(energy)
        axs[0].set_title('Total Energy Over Time')
        axs[0].set_xlabel('Time Steps')
        axs[0].set_ylabel('Energy (Joules)')
        axs[1].plot(np.linalg.norm(angular_momentum, axis=1))
        axs[1].set_title('Total Angular Momentum Magnitude Over Time')
        axs[1].set_xlabel('Time Steps')
        axs[1].set_ylabel(r'Angular Momentum $(kg m^2/s)$')
        plt.tight_layout()
        plt.show()