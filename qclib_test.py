import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator
from qclib.state_preparation import UCGInitialize
from qiskit.quantum_info import state_fidelity


# Generate 3-qubit random input state vector
n = 3
rnd = np.random.RandomState(42)
input_vector = rnd.rand(2 ** n) + rnd.rand(2 ** n) * 1j
input_vector = input_vector / np.linalg.norm(input_vector)

# Build a quantum circuit to initialize the input vector
circuit = UCGInitialize(input_vector).definition

# Construct an ideal simulator
backend = AerSimulator()

# Tests whether the produced state vector is equal to the input vector.
t_circuit = transpile(circuit, backend)
t_circuit.save_statevector()
state_vector = backend.run(t_circuit).result().get_statevector()
fidelity = state_fidelity(state_vector,input_vector)

print('Equal:', np.allclose(state_vector, input_vector))
print(f'Fidelity: {fidelity}')
# Equal: True
print(input_vector)
print(state_vector)