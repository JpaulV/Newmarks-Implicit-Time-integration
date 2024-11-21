# Newmark Method Simulation for 2-DOF System

This code simulates the dynamic response of a two-degree-of-freedom (2-DOF) system using the Newmark method (trapezoidal rule). The system consists of two masses connected by springs, and the displacements, velocities, and accelerations are computed over a specified time period.

## Overview
- **Author**: Jerry Paul Varghese, Adcom. Mat.N. 453553
- **Date**: 21/11/2024

The simulation uses the trapezoidal rule for time integration to compute the displacements, velocities, and accelerations of the system, and the results are plotted for visualization.

## Problem Description
- **System**: Two masses connected by springs with different stiffness values.
- **Stiffness and Mass Parameters**:
  - Mass at node 2: `m2 = 1`
  - Mass at node 3: `m3 = 1`
  - Stiffness between node 1 and node 2: `k1 = 10^7`
  - Stiffness between node 2 and node 3: `k2 = 1`
  - Prescribed frequency for displacement at node 1: `omega_p = 1.2`

## Newmark Parameters
- `beta = 0.25`
- `gamma = 0.5`

These parameters correspond to the trapezoidal rule for time integration.

## Implementation
- **Mass Matrix**: The mass matrix `M` is defined for nodes 2 and 3.
- **Stiffness Matrix**: The stiffness matrix `K` is defined for nodes 2 and 3.
- **Initial Conditions**: Initial displacements, velocities, and accelerations are set to zero.
- **Time Stepping**: The code calculates the response of the system over time using a prescribed displacement applied at node 1.
- **Plots**: The results are plotted in a single figure with four subplots showing:
  1. Displacements at nodes 2 and 3.
  2. Velocities at nodes 2 and 3.
  3. Acceleration at node 2.
  4. Acceleration at node 3.

## Usage
- **MATLAB Version**: This code is written in MATLAB and should be compatible with most recent versions.
- **Input Parameters**: You can modify the mass, stiffness, time step (`dt`), and simulation time (`total_time`) according to your problem requirements.

## Example Output
- The script generates plots for displacement, velocity, and acceleration for nodes 2 and 3 over the given simulation time. It provides insight into the dynamic response of the system, such as how the acceleration at node 3 flattens over time.

## Notes
- The prescribed displacement at node 1 is represented as an inverted sine wave (`-sin(omega_p * time(i))`).
- Make sure to adjust the time step size (`dt`) and the total simulation time (`total_time`) based on the desired accuracy and stability of the solution.

## License
This code is released under the MIT License.

Feel free to modify the code as needed for your project or research work.

## Contact
For any queries or suggestions, please contact Jerry Paul Varghese at `email@example.com`.
