from parser import extract_state_from_file
from isa import isa_prepare

# First argument is path to file, second argument is number of qubits
# Quantum state files are encoded in plain text. An [n]-qubit state is encoded
# by a file with 2^n lines, where each line has the format
#    [prefix] [amplitude] [phase]
# Such that the complex amplitude at basis state [prefix] is 
# [amplitude * e^(i * phase)]. Requires the [prefix]es to start from 0 and
# end at 2^n - 1, see "state.txt" as an example
state = extract_state_from_file("State.txt", 5)
gate_sequence = isa_prepare(state, target_fidelity=0.95)
CX_count = 0
for gate in gate_sequence:
    # print(gate.to_string())
    if gate.gate_type == "cx" or gate.gate_type == "xc":
        CX_count += 1
# print(f"Number of Gates: {len(gate_sequence)}")
print(f"CX Count: {CX_count}")
#[gate_sequence] is the list of gates to apply to state to approximately
# recover the |0> state:

import qc

state_c = state
for gate in gate_sequence:
    state_c = qc.apply_gate(gate, state_c)

print(abs(state_c[0]) ** 2) #Should be greater than 0.95

#To get the gate sequence for preparing [state] starting from |0>,
# reverse and invert all the gates

import numpy as np

start = [0 for _ in range(32)]
start[0] += 1
for gate in reversed(gate_sequence):
    start = qc.apply_gate(gate.inverse(), start)

start = np.array(start)
print(abs(np.vdot(state, start)) ** 2) #Should be greater than 0.95

