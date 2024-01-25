from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.providers.fake_provider import FakeMelbourneV2
from qiskit.quantum_info import state_fidelity
from qiskit.test.mock import FakeMelbourneV2

import numpy as np
import os
import qc
from parser import extract_state_from_file
import matplotlib.pyplot as plt
from datetime import datetime
import cmath


backend = AerSimulator.from_backend(FakeMelbourneV2())

# Retrieve the coupling map from the backend
print("Coupling Map:")
print(backend.configuration().coupling_map)
aer_sim = AerSimulator()
melbourne = FakeMelbourneV2()
n=5
state_file = 'state.txt'

state = extract_state_from_file(state_file, n)

print(state)
amplitudes_ideal = [abs(c) for c in state]

if "isa_cpp" not in os.listdir():
      os.system("./build_isa")
output_file  = "fidelity_gate_outputs/output.txt"
os.system(f"./isa_cpp {output_file} {state_file} {n} 0.80")

gate_sequence = []
with open(output_file) as f:
    while True:
        line = f.readline()
        if len(line) == 0: break
        if len(line) < 3: continue
        [gate_type, a1, a2] = line.split()
        if gate_type == "rx": 
            target = a1
            angle = a2
            gate = qc.Gate.RX(int(a1), float(angle), n)
            gate_sequence.append(gate.inverse())
        elif gate_type == "ry":
            target = a1
            angle = a2
            gate = qc.Gate.RY(int(a1), float(angle), n)                    
            gate_sequence.append(gate.inverse())
        elif gate_type == "rz":
            target = a1
            angle = a2
            gate = qc.Gate.RZ(int(a1), float(angle), n)
            gate_sequence.append(gate.inverse())
        elif gate_type == "cx" or gate_type == "xc":
            gate_sequence.append(qc.Gate.CX(int(a1), int(a2),n))

gate_sequence.reverse()

CX_count = 0
for gate in gate_sequence:
    # print(gate.to_string())
    if gate.gate_type == "cx" or gate.gate_type == "xc":
        CX_count += 1
print(f"Number of Gates: {len(gate_sequence)}")
print(f"CX Count: {CX_count}")



state_c = state
for gate in reversed(gate_sequence):
    state_c = qc.apply_gate(gate.inverse(), state_c)
ideal_fid1 = abs(state_c[0]) ** 2
print(f"Verify fid greater than 0.95: {ideal_fid1}") 

#To get the gate sequence for preparing [state] starting from |0>,
# reverse and invert all the gates

import numpy as np

start = [0 for _ in range(32)]
start[0] += 1
for gate in gate_sequence:
    start = qc.apply_gate(gate, start)

start = np.array(start)
ideal_fid2 = abs(np.vdot(state, start)) ** 2
print(f"Verify fid greater than 0.95: {ideal_fid2}")

circ  = QuantumCircuit(n)
for gate in gate_sequence:
    if gate.gate_type == "rx": 
        target = gate.target
        angle = gate.angle
        circ.rx(angle,target)
    elif gate.gate_type == "ry":
        target = gate.target
        angle = gate.angle
        circ.ry(angle,target)
    elif gate.gate_type == "rz":
        target = gate.target
        angle = gate.angle
        circ.rz(angle,target)

    elif gate.gate_type == "cx":
        control = gate.target
        target = control + 1
        circ.cx(control,target)

    elif gate.gate_type == "xc":
        target = gate.target
        control = target + 1
        circ.cx(control,target)
circ.save_density_matrix()
state = np.array(state).reshape(-1, 1)
p_ideal = np.outer(state,state.conj())
p_noisy = melbourne.run(circ).result().data()['density_matrix']
noisy_fidelity = state_fidelity(p_ideal,p_noisy)
print(f"Noisy fid : {noisy_fidelity}")
print(p_noisy)
probabilities = np.diag(p_noisy)
print(sum(probabilities))

# Calculate the square root to obtain the complex amplitudes
noisy_amplitudes = np.sqrt(probabilities)



# Plot both amplitudes and noisy amplitudes on the same plot
plt.plot(amplitudes_ideal, label='Ideal', linestyle='-')
plt.plot(noisy_amplitudes, label='Noisy', linestyle='-')

plt.title('Amplitudes Comparison (Ideal vs Noisy)')
plt.xlabel('Basis')
plt.ylabel('Amplitude')
plt.legend()
plt.savefig(f'Plots/5_qubit_amplitudes_0.8.png')
plt.show()

