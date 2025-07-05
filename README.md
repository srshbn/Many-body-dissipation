# Many-body-dissipation
Transverse-Field Ising Model with Dissipation (QuTiP Simulation)
This script simulates the dissipative dynamics of a transverse-field Ising model using the Lindblad master equation. It’s designed to explore parameter regimes for quantum experiments or emulations on near-term quantum devices.

# Description
Models a 1D chain of qubits with nearest-neighbor Ising interactions and a transverse magnetic field.
Includes local dissipation on each qubit using Lindblad formalism.
Computes time evolution of ⟨σz⟩ for each qubit.
Sweeps over different interaction-to-field ratios (J/h) to analyze dynamics.
Visualizes:
Expectation values over time
Dynamics of a selected qubit under various J/h
Final ⟨σz⟩ averaged over all qubits vs J/h
Heatmap of qubit 0 ⟨σz⟩ vs time and J/h
# Requirements
Install dependencies using:

pip install numpy matplotlib pandas qutip
How to Run
Run the script directly:

python ising_dissipative_simulation.py
# Steps 
Simulate the system dynamics
Plot key observables and trends under varying J/h
Optionally save expectation values as CSV
# Customize
Set N to change the number of qubits.
Modify J, h, and gamma to explore different physical regimes.
Use different observables or initial states by editing the simulation logic.
This simulation helps explore steady-state behavior and transient dynamics of open quantum systems before deploying real hardware runs.
