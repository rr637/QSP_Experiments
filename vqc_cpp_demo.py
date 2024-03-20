import os
import qc
from parser import extract_state_from_file
import numpy as np
from datetime import datetime
import time
import csv

#Compile cpp code to executable, if it hasn't been compiled already.
# Or, if there's an issue with running the executable, try removing it and
# recompiling
num_layers = 71
fidelity = 0
params = {'protein': 'homo_sapien', 'num_protein': 20, 'num_qubits':10, 'Avg_CX': 0, 'Avg_Runtime(s)': 0}
n = 10
CX_vqc_counts = []
Runtimes_vqc = []
Fidelities = []
name_exp = datetime.now().strftime("%m_%d_%H_%M_%S")
directory_path = f'Results/VQC_Exp:_{name_exp}'
output_vqc_path = f'Results/VQC_Exp:_{name_exp}/Outputs_vqc'
os.makedirs(output_vqc_path, exist_ok=True)

for i in range(params['num_protein']):
    
    if "vqc_cpp" in os.listdir():
        os.remove("vqc_cpp")

    if "vqc_cpp" not in os.listdir():
        os.system("./build_vqc")

    state_file = f"./States/state{i}.txt"
    output_file = f"Results/VQC_Exp:_{name_exp}/Outputs_vqc/gates{i}.txt"
    start_time = time.time()

    os.system(f"./vqc_cpp {output_file} {state_file} {n} 0.95 {num_layers}")
    runtime = time.time() - start_time
    state = extract_state_from_file(f"{state_file}", n)
    gate_sequence = []
    with open(output_file) as f:
        while True:
            line = f.readline()
            if len(line) == 0: break
            if len(line) < 3: continue
            [gate_type, a1, a2] = line.split()
            if gate_type == "rx":
                gate_sequence.append(qc.Gate.RX(int(a1), float(a2), n))
            elif gate_type == "ry":
                gate_sequence.append(qc.Gate.RY(int(a1), float(a2), n))
            elif gate_type == "rz":
                gate_sequence.append(qc.Gate.RZ(int(a1), float(a2), n))
            elif gate_type == "cx" or gate_type == "xc":
                gate_sequence.append(qc.Gate.CX(int(a1), int(a2), n))
            else: raise

    CX_count = 0
    for gate in gate_sequence:
        # print(gate.to_string())
        if gate.gate_type == "cx" or gate.gate_type == "xc":
            CX_count += 1
    print(f"Number of Layers: {num_layers}")
    print(f"Number of Gates: {len(gate_sequence)}")
    print(f"CX Count: {CX_count}")
    CX_vqc_counts.append(CX_count)
    Runtimes_vqc.append(runtime)

    state_c = state
    for gate in reversed(gate_sequence):
        state_c = qc.apply_gate(gate.inverse(), state_c)
    fidelity = abs(state_c[0]) ** 2
    print(f"Fidelity: {fidelity}") #Should be greater than 0.95
    Fidelities.append(fidelity)
    #To get the gate sequence for preparing [state] starting from |0>,
    # reverse and invert all the gates

    start = [0 for _ in range(2**n)]
    start[0] += 1
    for gate in gate_sequence:
        start = qc.apply_gate(gate, start)

    start = np.array(start)
    print(abs(np.vdot(state, start)) ** 2) #Should be greater than 0.95

avg_cx = sum(CX_vqc_counts) / len(CX_vqc_counts)
avg_runtimes = sum(Runtimes_vqc) / len(Runtimes_vqc)
avg_fidelities = sum(Fidelities) / len(Fidelities)
csv_file_path = f'{directory_path}/experiment_results.csv'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['CX_counts'] + CX_vqc_counts)
    writer.writerow(['Avg_CX', avg_cx])
    writer.writerow(['Runtimes'] + Runtimes_vqc)
    writer.writerow(['Avg_Runtime', avg_runtimes])
    writer.writerow(['Fidelities'] + Fidelities)
    writer.writerow(['Avg_Fidelities', avg_fidelities])



# while fidelity < 0.95: 
#     if "vqc_cpp" in os.listdir():
#         os.remove("vqc_cpp")

#     if "vqc_cpp" not in os.listdir():
#         os.system("./build_vqc")

#     n = 10
#     state_file = "States/state0.txt"
#     #Same command as running the ISA executable. To change the number of VQC layers,
#     # go into cpp/vqc_main.cpp and modify the variable "num_layers" under main()
#     os.system(f"./vqc_cpp output4.txt {state_file} {n} 0.95 {num_layers}")
#     state = extract_state_from_file(f"{state_file}", n)

#     #Read the gates from the output file, same as for ISA
#     gate_sequence = []
#     with open("output4.txt") as f:
#         while True:
#             line = f.readline()
#             if len(line) == 0: break
#             if len(line) < 3: continue
#             [gate_type, a1, a2] = line.split()
#             if gate_type == "rx":
#                 gate_sequence.append(qc.Gate.RX(int(a1), float(a2), n))
#             elif gate_type == "ry":
#                 gate_sequence.append(qc.Gate.RY(int(a1), float(a2), n))
#             elif gate_type == "rz":
#                 gate_sequence.append(qc.Gate.RZ(int(a1), float(a2), n))
#             elif gate_type == "cx" or gate_type == "xc":
#                 gate_sequence.append(qc.Gate.CX(int(a1), int(a2), n))
#             else: raise

#     CX_count = 0
#     for gate in gate_sequence:
#         # print(gate.to_string())
#         if gate.gate_type == "cx" or gate.gate_type == "xc":
#             CX_count += 1
#     print(f"Number of Layers: {num_layers}")
#     print(f"Number of Gates: {len(gate_sequence)}")
#     print(f"CX Count: {CX_count}")


#     state_c = state
#     for gate in reversed(gate_sequence):
#         state_c = qc.apply_gate(gate.inverse(), state_c)
#     fidelity = abs(state_c[0]) ** 2
#     print(f"Fidelity: {fidelity}") #Should be greater than 0.95

#     #To get the gate sequence for preparing [state] starting from |0>,
#     # reverse and invert all the gates

#     start = [0 for _ in range(2**n)]
#     start[0] += 1
#     for gate in gate_sequence:
#         start = qc.apply_gate(gate, start)

#     start = np.array(start)
#     print(abs(np.vdot(state, start)) ** 2) #Should be greater than 0.95

#     num_layers += 2

