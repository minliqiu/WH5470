# WH5470

Final Project for ASTR 5470: Wisdom-Holman N-body Symplectic Integrator

This project simulates the gravitational interactions between multiple bodies in space using the Wisdom-Holman symplectic integrator, an efficient method for long-term numerical simulations of planetary systems. It's designed to conserve the system's Hamiltonian over time, making it practical for studying orbital dynamics.

---

## Features

- **Symplectic Integration**: Uses the Wisdom-Holman method to ensure better long-term energy conservation.
- **3D Visualization**: Includes functionality to visualize the trajectories of celestial bodies in three dimensions.
- **Energy and Angular Momentum Conservation Check**: Computes and plots total energy and angular momentum to verify their conservation over the simulation period.
- **Close Encounter Detection**: Integrations are paused when a hill radius close encounter happens.
- **General Relativity Modification**: Adds optional gr module to modify general relativity effects.
- **Long Term Evolution Stability Check**: Wisdom-Holman integrators are best used for long-term evolution and stability check.

  
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
git clone https://github.com/ASTRX470/final-project-minliqiu.git
```
### Usage

Navigate to the project directory and run the script:

```bash
cd n-body-simulation
python wh5470.py
```

This will start a simulation with predefined initial conditions and open a plot showing the trajectories of the bodies.

## Configuration

To change the simulation parameters (e.g., number of bodies, time step), modify the parameters directly in `wh5470.py`.

## Authors

- **Minli Qiu** - *Initial work* - [minliqiu](https://github.com/minliqiu/wh5470)

---
