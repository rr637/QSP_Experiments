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
params = {'protein': 'homo_sapien', 'num_protein': 20, 'num_qubits':10, 'Avg_CX': 0, 'Avg_Runtime(s)': 0}
n = 10
CX_isa_counts = []
Runtimes_isa = []
Fidelities = []
name_exp = datetime.now().strftime("%m_%d_%H_%M_%S")
directory_path = f'Results/ISA_Exp:_{name_exp}'
output_vqc_path = f'Results/ISA_Exp:_{name_exp}/Outputs_isa'
os.makedirs(output_vqc_path, exist_ok=True)

for i in range(params['num_protein']):
    
    if "isa_cpp" in os.listdir():
        os.remove("isa_cpp")

    if "isa_cpp" not in os.listdir():
        os.system("./build_isa")

    state_file = f"./States/state{i}.txt"
    output_file = f"Results/ISA_Exp:_{name_exp}/Outputs_isa/gates{i}.txt"
    start_time = time.time()

    os.system(f"./isa_cpp {output_file} {state_file} {n} 0.95")
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
    print(f"Number of Gates: {len(gate_sequence)}")
    print(f"CX Count: {CX_count}")
    CX_isa_counts.append(CX_count)
    Runtimes_isa.append(runtime)

    state_c = state
    for gate in gate_sequence:
        state_c = qc.apply_gate(gate, state_c)
    fidelity = abs(state_c[0]) ** 2
    print(f"Fidelity: {fidelity}") #Should be greater than 0.95
    Fidelities.append(fidelity)

    #To get the gate sequence for preparing [state] starting from |0>,
    # reverse and invert all the gates

    start = [0 for _ in range(2**n)]
    start[0] += 1
    for gate in reversed(gate_sequence):
        start = qc.apply_gate(gate.inverse(), start)

    start = np.array(start)
    print(abs(np.vdot(state, start)) ** 2) #Should be greater than 0.95

avg_cx = sum(CX_isa_counts) / len(CX_isa_counts)
avg_runtimes = sum(Runtimes_isa) / len(Runtimes_isa)
avg_fidelities = sum(Fidelities) / len(Fidelities)
csv_file_path = f'{directory_path}/experiment_results.csv'
with open(csv_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['CX_counts'] + CX_isa_counts)
    writer.writerow(['Avg_CX', avg_cx])
    writer.writerow(['Runtimes'] + Runtimes_isa)
    writer.writerow(['Avg_Runtime', avg_runtimes])
    writer.writerow(['Fidelities'] + Fidelities)
    writer.writerow(['Avg_Fidelities', avg_fidelities])

import os
import qc

#Compile cpp code to executable, if it hasn't been compiled already
if "isa_cpp" not in os.listdir():
    os.system("./build_isa")

#In the below command,
#  [output.txt] is the output file (can be named anything)
#  [state.txt] is the input file
#  [5] is the number of qubits
#  [0.95] is the target fidelity. The target fidelity is an optional argument
#    with 0.95 as default value
os.system("./isa_cpp output.txt state.txt 5 0.95")

#Now read the gates from the output file
gate_sequence = []
with open("output.txt") as f:
    while True:
        line = f.readline()
        if len(line) == 0: break
        if len(line) < 3: continue
        [gate_type, a1, a2] = line.split()
        if gate_type == "rx": 
            gate_sequence.append(qc.Gate.RX(int(a1), float(a2), 5))
        elif gate_type == "ry":
            gate_sequence.append(qc.Gate.RY(int(a1), float(a2), 5))
        elif gate_type == "rz":
            gate_sequence.append(qc.Gate.RZ(int(a1), float(a2), 5))
        elif gate_type == "cx" or gate_type == "xc":
            gate_sequence.append(qc.Gate.CX(int(a1), int(a2), 5))
        else: raise

#[gate_sequence] is the same as 
# gate_sequence = isa_prepare(state, target_fidelity=0.95)
# from py_demo.py
# Getting the gate sequence using the C++ executable is a bit more convoluted
# but the C++ code is much much faster than the Python version.







# #Compile cpp code to executable, if it hasn't been compiled already
# if "isa_cpp" not in os.listdir():
#     os.system("./build_isa")
# n = 10
# #In the below command,
# #  [output.txt] is the output file (can be named anything)
# #  [state.txt] is the input file
# #  [5] is the number of qubits
# #  [0.95] is the target fidelity. The target fidelity is an optional argument
# #    with 0.95 as default value
# state_file = 'States/state0.txt'
# os.system(f"./isa_cpp output1.txt {state_file} {n} 0.95")

# #Now read the gates from the output file
# gate_sequence = []
# with open("output3.txt") as f:
#     while True:
#         line = f.readline()
#         if len(line) == 0: break
#         if len(line) < 3: continue
#         [gate_type, a1, a2] = line.split()
#         if gate_type == "rx": 
#             gate_sequence.append(qc.Gate.RX(int(a1), float(a2), n ))
#         elif gate_type == "ry":
#             gate_sequence.append(qc.Gate.RY(int(a1), float(a2), n))
#         elif gate_type == "rz":
#             gate_sequence.append(qc.Gate.RZ(int(a1), float(a2), n))
#         elif gate_type == "cx" or gate_type == "xc":
#             gate_sequence.append(qc.Gate.CX(int(a1), int(a2), n))
#         else: raise

# #[gate_sequence] is the same as 
# # gate_sequence = isa_prepare(state, target_fidelity=0.95)
# # from py_demo.py
# # Getting the gate sequence using the C++ executable is a bit more convoluted
# # but the C++ code is much much faster than the Python version.

# CX_count = 0
# for gate in gate_sequence:
#     # print(gate.to_string())
#     if gate.gate_type == "cx" or gate.gate_type == "xc":
#         CX_count += 1
# print(f"Number of Gates : {len(gate_sequence)}")
# print(f"CX Count: {CX_count}")
# print(state_file)