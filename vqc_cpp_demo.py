import os
import qc
from parser import extract_state_from_file
import numpy as np
#Compile cpp code to executable, if it hasn't been compiled already.
# Or, if there's an issue with running the executable, try removing it and
# recompiling
num_layers = 6
fidelity = 0
while fidelity < 0.95: 
    if "vqc_cpp" in os.listdir():
        os.remove("vqc_cpp")

    if "vqc_cpp" not in os.listdir():
        os.system("./build_vqc")

    n = 5
    state_file = "state.txt"
    #Same command as running the ISA executable. To change the number of VQC layers,
    # go into cpp/vqc_main.cpp and modify the variable "num_layers" under main()
    os.system(f"./vqc_cpp output4.txt {state_file} {n} 0.95 {num_layers}")
    state = extract_state_from_file(f"{state_file}", n)

    #Read the gates from the output file, same as for ISA
    gate_sequence = []
    with open("output4.txt") as f:
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


    state_c = state
    for gate in reversed(gate_sequence):
        state_c = qc.apply_gate(gate.inverse(), state_c)
    fidelity = abs(state_c[0]) ** 2
    print(f"Fidelity: {fidelity}") #Should be greater than 0.95

    #To get the gate sequence for preparing [state] starting from |0>,
    # reverse and invert all the gates

    start = [0 for _ in range(2**n)]
    start[0] += 1
    for gate in gate_sequence:
        start = qc.apply_gate(gate, start)

    start = np.array(start)
    print(abs(np.vdot(state, start)) ** 2) #Should be greater than 0.95

    num_layers += 2

