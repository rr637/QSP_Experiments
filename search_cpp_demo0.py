import os
import qc
import time

#Compile cpp code to executable, if it hasn't been compiled already.
# Or, if there's an issue with running the executable, try removing it and
# recompiling
if "search_cpp" not in os.listdir():
    os.system("./build_search")

n = 10

start_time = time.time()

#Same command as running the ISA executable. To change the number of VQC layers,
# go into cpp/vqc_main.cpp and modify the variable "num_layers" under main()
state_file = 'States/state0.txt'
os.system(f"./search_cpp output.txt {state_file} {n} 0.95")

#Read the gates from the output file, same as for ISA
gate_sequence = []
with open("output.txt") as f:
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

end_time = time.time()
runtime = end_time - start_time
CX_count = 0
for gate in gate_sequence:
    # print(gate.to_string())
    if gate.gate_type == "cx" or gate.gate_type == "xc":
        CX_count += 1

        
print(f"Number of Gates: {len(gate_sequence)}")
print(f"CX Count: {CX_count}")
print(f"Script runtime: {runtime} seconds")
print(state_file)