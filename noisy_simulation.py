from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeProvider, FakeMelbourneV2
from qiskit.quantum_info import state_fidelity
import numpy as np
import os
import qc
from parser import extract_state_from_file
import matplotlib.pyplot as plt
from datetime import datetime
import statistics

backend = FakeMelbourneV2()

# Retrieve the coupling map from the backend

aer_sim = AerSimulator()
melbourne = FakeMelbourneV2()
n=5


def run_exp(state):
  backend = AerSimulator.from_backend(FakeMelbourneV2())

# Retrieve the coupling map from the backend
  print("Coupling Map:")
  print(backend.configuration().coupling_map)
  aer_sim = AerSimulator()
  melbourne = FakeMelbourneV2()
  n=5
  state_file = "5_qubit_state.txt"
  with open(f"{state_file}", "w") as f:
        for prefix, amp, phase in state:
            f.write(f"{prefix} {amp} {phase}\n")
  state = extract_state_from_file(state_file,n)
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
    target_fids[str(CX_count)] = ideal_fid2

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
  print(f"CX-Count: idealfid - {target_fids}")
  print(f"CX-Count: noisyfid - {noisy_fids}")
  return (target_fids,noisy_fids, noisy_vs_target)

import random
import cmath

def rand_state(n):
    state = np.array([random.gauss(0, 1) for _ in range(1 << n)])
    norm = state.dot(state) ** 0.5
    return state / norm

def rand_complex_state(n):
    amps = rand_state(n)
    phases = np.array([cmath.exp(2.0j * cmath.pi * random.random()) for _ in range(1 << n)])
    return amps * phases

def complex_to_amp_phase(c):
    amp = abs(c)
    phase = cmath.phase(c)
    return amp, phase

def format_number(number):
    # Format to 3 significant figures
    return '{:.6g}'.format(number)


target_fids_per_cx = {str(i): [] for i in range(0, 51)}

noisy_fids_per_cx =  {str(i): [] for i in range(0, 51)}
noisy_vs_target_points = []
n_qubits = 5
n_states = 100
def add_wo_dup(list,e,threshold):
    for item in list:
        if abs(item - e) <= threshold:
            return False  # Duplicate found within threshold, don't add
    list.append(e)  # No duplicates found, add the new element
    return True
threshold = 1e-4
for i in range(n_states):
    state = rand_complex_state(n_qubits)

    # Convert complex numbers to amplitude and phase
    state = [(format(i, f'0{n_qubits}b'), format_number(amp), format_number(phase)) 
                    for i, (amp, phase) in enumerate(map(complex_to_amp_phase, state))]

    # Write state to a text file
    with open("output.txt", "w") as f:
        for prefix, amp, phase in state:
            f.write(f"{prefix} {amp} {phase}\n")


    target_fids, noisy_fids, noisy_vs_target = run_exp(state)
    for cx,fid in target_fids.items():
        add_wo_dup(target_fids_per_cx[cx],fid,threshold)
    for cx,fid in noisy_fids.items():
        add_wo_dup(noisy_fids_per_cx[cx],fid,threshold)
    for target,noisy in noisy_vs_target.items():
        noisy_vs_target_points.append((float(target),noisy))


# Remove elements with empty lists from target_fids_per_cx
target_fids_per_cx = {key: value for key, value in target_fids_per_cx.items() if value}

# Remove elements with empty lists from noisy_fids_per_cx
noisy_fids_per_cx = {key: value for key, value in noisy_fids_per_cx.items() if value}

print(f"tfdict: {target_fids_per_cx}")
print(f"nfdict: {noisy_fids_per_cx}")
print(f"ntpoints: {noisy_vs_target_points}")

def sd(list):
    if len(list)<2:
        return  0
    else:
        return statistics.stdev(list)

target_avg_err  = {cx:(sum(list)/len(list),sd(list)) for cx,list in target_fids_per_cx.items()}
noisy_avg_err  = {cx:(sum(list)/len(list), sd(list)) for cx,list in noisy_fids_per_cx.items()}
print(f"target_means,error: {target_avg_err}")
print(f"noisy_means,error: {noisy_avg_err}  ")

x_values_target = [int(key) for key in target_avg_err.keys()]
y_values_target, std_target = zip(*target_avg_err.values())

x_values_noisy = [int(key) for key in noisy_avg_err.keys()]
y_values_noisy, std_noisy = zip(*noisy_avg_err.values())

x_points,y_points = zip(*noisy_vs_target_points)


# Create a 1x2 grid
fig, axs = plt.subplots(2, 2, figsize=(12, 6))
y_ticks = np.arange(0, 1.1, 0.2)

# Plot for noisy_avg_err with error bars
axs[0,0].bar(x_values_noisy, y_values_noisy, yerr=std_noisy, color='orange', alpha=0.7)
axs[0,0].set_xlabel('CX Count',fontsize = 13)
axs[0,0].set_ylabel('ISA - Noisy Fidelity ',fontsize = 13)
axs[0,0].set_yticks(y_ticks)
axs[0,0].set_ylim(0, 1.05)

# Plot for target_avg_err with error bars
axs[0,1].bar(x_values_target, y_values_target, yerr=std_target, color='blue', alpha=0.7)
axs[0,1].set_xlabel('CX Count',fontsize = 13)
axs[0,1].set_ylabel('ISA - Theoretical Fidelity',fontsize = 13)
axs[0,1].set_yticks(y_ticks)
axs[0,1].set_ylim(0, 1.05)


axs[1, 0].scatter(x_points, y_points, color='red', s=5)  # Adjust the value of 's' as needed for smaller points
axs[1,0].set_xlabel('ISA - Theoretical Fidelty',fontsize = 13)
axs[1,0].set_ylabel('ISA - Noisy Fidelity',fontsize = 13)
axs[1,0].set_ylim(0.3, 0.71)
axs[1,0].set_xlim(0.3,1.05)
axs[1,0].set_yticks(np.arange(0.3,0.71,0.1))

plt.tight_layout()



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









# Assuming you have defined amplitudes_ideal and noisy_amplitudes for the bottom-right subplot

# Plot for amplitudes comparison in the fourth subplot (bottom-right)
tick_positions = list(range(0, 33, 4))  # Show tick marks every 4 ticks, including 0 and 32

# Define the basis states for ticks
basis_states = ['' for _ in range(9)]  # Initialize all ticks with empty strings
basis_states[0] = r'$|\text{00000}\rangle$'  # Set label at x=0
basis_states[8] = r'$|\text{11111}\rangle$'  # Set label at x=32
basis_states[2] = r'$\dots$'  # Set ellipsis at x=8
basis_states[4] = r'$\dots$'  # Set ellipsis at x=16
basis_states[6] = r'$\dots$'  # Set ellipsis at x=24

# Set ticks and labels
axs[1, 1].set_xticks(tick_positions)
axs[1, 1].set_xticklabels(basis_states, fontsize=10)  # Set font size

# Plot ideal and noisy amplitudes
axs[1, 1].plot(amplitudes_ideal, label='Target Amplitudes', linestyle='-', marker='o', markersize=2)
axs[1, 1].plot(noisy_amplitudes, label='Noisy Amplitudes', linestyle='-', marker='o', markersize=2)


# Set labels and legend
axs[1, 1].set_xlabel('Basis States', fontsize = 13)
axs[1, 1].set_ylabel('$|\mathrm{Amplitude}|^2$', fontsize=13)
axs[1, 1].legend()


# Add ellipsis at tick 16

# Adjust layout for better appearance
plt.tight_layout()
current_time = datetime.now().strftime('%Y%m%d_%H%M%S')
plt.savefig(f'Plots/Final_wo_dup{current_time}.png')
plt.show()




# #{'4': 0.4156102611680524, '10': 0.5617238763523772, '14': 0.5741604922403004, '19': 0.5638220750729754, '24': 0.5218784637169586, '30': 0.4934750850822063}
# #{'4': 0.4156102611680524, '10': 0.5617238763523772, '14': 0.5741604922403004, '19': 0.5638220750729754, '24': 0.5218784637169586, '30': 0.4934750850822063}