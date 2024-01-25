from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.providers.fake_provider import FakeMelbourneV2
from qiskit.quantum_info import state_fidelity
import numpy as np
import os
import qc
from parser import extract_state_from_file
import matplotlib.pyplot as plt
from datetime import datetime

backend = AerSimulator.from_backend(FakeMelbourneV2())

# Retrieve the coupling map from the backend
print("Coupling Map:")
print(backend.configuration().coupling_map)
aer_sim = AerSimulator()
melbourne = FakeMelbourneV2()
n=5
state = extract_state_from_file("state.txt", n)



def run_exp():
  backend = AerSimulator.from_backend(FakeMelbourneV2())

# Retrieve the coupling map from the backend
  print("Coupling Map:")
  print(backend.configuration().coupling_map)
  aer_sim = AerSimulator()
  melbourne = FakeMelbourneV2()
  n=5
  state = extract_state_from_file("state.txt", n)

  #Compile cpp code to executable, if it hasn't been compiled already
  if "isa_cpp" not in os.listdir():
      os.system("./build_isa")

  target_fids_list = ['0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8','0.85','0.90','0.95','0.98']
  CX_counts  = {}
  noisy_fids = {}
  target_fids = {}
  noisy_vs_target = {}
  #In the below command,
  #  [output.txt] is the output file (can be named anything)
  #  [state.txt] is the input file
  #  [5] is the number of qubits
  #  [0.95] is the target fidelity. The target fidelity is an optional argument
  #    with 0.95 as default value
  state_file = 'state.txt'
  for s in target_fids_list:
    os.system(f"./isa_cpp fidelity_gate_outputs/output_{s}.txt {state_file} {n} {s}")

    #Now read the gates from the output file
    gate_sequence = []
    with open(f"fidelity_gate_outputs/output_{s}.txt") as f:
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
      
    #[gate_sequence] is the same as 
    # gate_sequence = isa_prepare(state, target_fidelity=0.95)
    # from py_demo.py
    # Getting the gate sequence using the C++ executable is a bit more convoluted
    # but the C++ code is much much faster than the Python version.

    CX_count = 0
    for gate in gate_sequence:
        # print(gate.to_string())
        if gate.gate_type == "cx" or gate.gate_type == "xc":
            CX_count += 1
    CX_counts[s]  = CX_count
    target_fids[str(CX_count)] = float(s)
    state_c = state
    for gate in reversed(gate_sequence):
        state_c = qc.apply_gate(gate.inverse(), state_c)
    ideal_fid1 = abs(state_c[0]) ** 2
    print(f"Verify fid greater than {s}: {ideal_fid1}") 

    #To get the gate sequence for preparing [state] starting from |0>,
    # reverse and invert all the gates

    import numpy as np

    start = [0 for _ in range(32)]
    start[0] += 1
    for gate in gate_sequence:
        start = qc.apply_gate(gate, start)

    start = np.array(start)
    ideal_fid2 = abs(np.vdot(state, start)) ** 2
    print(f"Verify fid greater than {s}: {ideal_fid2}") 
    CX_count = 0
    for gate in gate_sequence:
        # print(gate.to_string())
        if gate.gate_type == "cx" or gate.gate_type == "xc":
            CX_count += 1
    CX_counts[ideal_fid2]  = CX_count
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
    fidelity = state_fidelity(p_ideal,p_noisy)
    noisy_fids[str(CX_count)] = fidelity
    noisy_vs_target[str(ideal_fid2)] = fidelity
  print(target_fids)
  print(noisy_fids)
  return (target_fids,noisy_fids, noisy_vs_target)


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








# Assuming target_fids, noisy_fids, and noisy_vs_target are dictionaries
target_fids, noisy_fids, noisy_vs_target = run_exp()

# Assuming you have already defined target_fids, noisy_fids, and noisy_vs_target

# Create a 2x2 grid
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Plot for noisy_fids as a bar graph in the first subplot (top-left)
x_values_noisy = [float(key) for key in noisy_fids.keys()]
y_values_noisy = list(noisy_fids.values())
axs[0, 0].bar(x_values_noisy, y_values_noisy, color='orange')
axs[0, 0].set_title('a) Noisy Fidelity')
axs[0, 0].set_xlabel('CX Count')
axs[0, 0].set_ylabel('Noisy Fidelity')

# Plot for target_fids as a bar graph in the second subplot (top-right)
x_values_target = [float(key) for key in target_fids.keys()]
y_values_target = list(target_fids.values())
axs[0, 1].bar(x_values_target, y_values_target, color='blue')
axs[0, 1].set_title('b) Ideal Fidelity')
axs[0, 1].set_xlabel('CX Count')
axs[0, 1].set_ylabel('Ideal Fidelity')

# Plot for noisy_vs_target as a line graph in the third subplot (bottom-left)
x_values_nvt = [float(key) for key in noisy_vs_target.keys()]
y_values_nvt = list(noisy_vs_target.values())
axs[1, 0].plot(x_values_nvt, y_values_nvt, color='green', marker='o')
axs[1, 0].set_title('c) Noisy vs Ideal Fidelity')
axs[1, 0].set_xlabel('Ideal Fidelity')
axs[1, 0].set_ylabel('Noisy Fidelity')

# Assuming you have defined amplitudes_ideal and noisy_amplitudes for the bottom-right subplot
# Plot for amplitudes comparison in the fourth subplot (bottom-right)
axs[1, 1].plot(amplitudes_ideal, label='Ideal', linestyle='-')
axs[1, 1].plot(noisy_amplitudes, label='Noisy', linestyle='-')
axs[1, 1].set_title('d) Amplitudes Comparison (Ideal vs Noisy)')
axs[1, 1].set_xlabel('Basis')
axs[1, 1].set_ylabel('Amplitude')
axs[1, 1].legend()

# Adjust layout for better appearance
plt.tight_layout()

# Save the figure
current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
plt.savefig(f'Plots/2x2_subplot_grid_{current_time}.png')

# Show the plot
plt.show()
