# WH5470

Final Project for ASTR 5470: Wisdom-Holman N-body Symplectic Integrator

This project simulates the gravitational interactions between multiple bodies in space using the Wisdom-Holman symplectic integrator, an efficient method for long-term numerical simulations of planetary systems. It's designed to conserve the system's Hamiltonian over time, making it practical for studying orbital dynamics.

---

## Features

- **Symplectic Integration**: Uses the Wisdom-Holman method to ensure better long-term energy conservation.
- **3D Visualization**: Includes functionality to visualize the trajectories of celestial bodies in three dimensions.
- **Energy and Angular Momentum Conservation Check**: Computes and plots total energy and angular momentum to verify their conservation over the simulation period.
- **Close Encounter Detection**:
- **General Relativity Modification**:
- **Long Term Evolution Stability Check**:

  
## Getting Started

### Prerequisites

- Python 3.x
- NumPy
- Matplotlib
- mpl_toolkits for 3D plotting
- rebound (optional, only required for Test 3)

### Installation

Clone this repository to your local machine using:

```bash
git clone https://github.com/your_username/n-body-simulation.git
```

### Usage

Navigate to the project directory and run the script:

```bash
cd n-body-simulation
python n_body_simulator.py
```

This will start a simulation with predefined initial conditions and open a plot showing the trajectories of the bodies.

## Configuration

To change the simulation parameters (e.g., number of bodies, time step), modify the `config.json` file or adjust the parameters directly in `n_body_simulator.py`.

## Contributing

Contributions are welcome! For major changes, please open an issue first to discuss what you would like to change. Please ensure to update tests as appropriate.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

- **Your Name** - *Initial work* - [YourUsername](https://github.com/YourUsername)

## Acknowledgments

- Hat tip to anyone whose code was used
- Inspiration
- etc

---
