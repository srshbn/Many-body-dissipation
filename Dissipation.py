import numpy as np
from qutip import *
import pandas as pd
import matplotlib.pyplot as plt

"""
This code simulates transverse-field Ising model with dissipative dynamics to optimize the parameters before running on a quantum computer: 

Note: A separate simulation handles **non-equilibrium** and **auxiliary-driven dynamics**, which is included in the package. 

It models a chain of N interacting qubits with nearest-neighbor Ising interactions and a transverse magnetic field,
with local dissipation (modeled in the Lindblad master equation)

Hamiltonian is given as the following:

  H = -J sum(sigmaxk sigmax(k+1) - h sum sigmazk
    - J: Coupling strength between neighbors (Ising term)
    - h: Transverse field strength

Local dissipation is modeled by collapse operators using gamma as dissipation rate and annihilation operator

The code performs the following:

1. Simulation of time evolution using the Lindblad master equation
2. It computes the evolution of ⟨sigmaz⟩ for each qubit over time
3. The dynamics of Qubit 1 under different ratios of interaction strength to field strength
4. The final average ⟨sigmaz⟩ across all qubits as a function of J/h
5. A heatmap of ⟨sigmaz⟩ for qubit 0 over time and varying J/h ratios

"""

# Parameters
N = 4         # Number of qubits
#Note: the runtime depends on N as Hilbert space tensors are 2**N and every operation acts on this entire space
J = 0.7         # Interaction strength
h = 1          # Transverse field strength
gamma = 0.5      # Dissipation rate
tlist = np.linspace(0, 10, 500)  # Time list

def simulate_dissipation(N, J, h, gamma):
    """
        Simulates the dissipative in transverse-field Ising model for N qubits.

        Parameters:
        - N: Number of qubits
        - J: Coupling strength between neighbors (Ising term)
        - h: Transverse field strength
        - gamma: Dissipation rate
    """

    # Pauli matrices
    sx = sigmax()
    sz = sigmaz()
    si = qeye(2) # identity operator [[1,0],[0,1]]
    #  Hamiltonian of transverse-field Ising model
    H = sum(-J * tensor([si] * i + [sx, sx] + [si] * (N - i - 2)) for i in range(N - 1))
    H += sum(-h * tensor([si] * i + [sz] + [si] * (N - i - 1)) for i in range(N))
    # Collapse operators for dissipation
    c_ops = [np.sqrt(gamma) * tensor([si] * i + [destroy(2)] + [si] * (N - i - 1)) for i in range(N)]
    # Initial state of qubits: all spins down
    psi0 = tensor([basis(2, 0)] * N)
    # Solve the master equation using mesolve and applying the time evolution
    result = mesolve(H, psi0, tlist, c_ops=c_ops, e_ops=[])
    # Expectation values after measurment
    expect_sz = np.array([expect(tensor([si] * i + [sz] + [si] * (N - i - 1)), result.states) for i in range(N)])
    # expect_sz shape: (N qubits, time steps)
    # Transpose it to DataFrame so each row = time, each column = qubit
    df = pd.DataFrame(expect_sz.T, columns=[f'Qubit {i}' for i in range(N)])
    # Add the time list as a new column
    df['Time'] = tlist
    df.set_index('Time', inplace=True) # Setting the time as the index of the DataFrame
    return df

# Calling the function for the desired parameters
df = simulate_dissipation(N=4, J=0.7, h=1, gamma=0.5)

#-------------------------------------------------------------------------------------------------------------

# Plot the results
df.plot(figsize=(10, 6), title='Expectation Value sigma_z for each qubit as a function of time')
plt.xlabel('Time')
plt.ylabel('Expectation value of sigma_z')
plt.legend()
plt.title('Dissipative dynamics of the transverse-field Ising model')
# plot qubit 1: ratios between J and h
ratios = [0.1, 0.5, 1.0, 1.5, 2.0] # J/h
colors = ['blue', 'green', 'orange', 'red', 'purple']
plt.figure()
for r, c in zip(ratios, colors):
    df = simulate_dissipation(N=4, J=r, h=1, gamma=0.5)
    plt.plot(df.index, df['Qubit 1'], label=f'J/h = {r}', color = c)
plt.xlabel('Time')
plt.ylabel('Expectation value of sigma_z')
plt.legend()
# plot sigma_z vs J/h at final time
final_vals = []
ratios = np.linspace(0.1, 2.0, 20)  # smoother sweep
for r in ratios:
    df =simulate_dissipation(N=4, J=r, h=1.0, gamma=0.5)
    final_sz = df.iloc[-1].mean()  # mean ⟨σz⟩ across all qubits at final time
    final_vals.append(final_sz)
plt.figure(figsize=(8, 5))
plt.plot(ratios, final_vals, marker='o')
plt.title("Final ⟨sigma_z⟩ (averaged over qubits) vs J/h")
plt.xlabel("J/h")
plt.ylabel("Final ⟨sigma_z⟩ at t = 10")
plt.grid(True)
#plot the heatmap
heatmap_data = []
for r in ratios:
    J = r*h
    sz_vals = simulate_dissipation(N, J, h, gamma)
    heatmap_data.append(sz_vals['Qubit 0'])
heatmap_data = np.array(heatmap_data)
plt.figure(figsize=(8, 5))
plt.imshow(
    heatmap_data,
    aspect='auto',
    cmap='coolwarm',
    extent=[tlist[0], tlist[-1], ratios[-1], ratios[0]],
    vmin=0,
    vmax=0.8
)
plt.colorbar(label='⟨σz⟩ on Qubit 0')
plt.title("Heatmap of ⟨σz⟩ on Qubit 0 vs Time and J/h")
plt.xlabel("Time")
plt.ylabel("J/h")
plt.tight_layout()

plt.show()

#df.to_csv("qubit_expectations.csv")
