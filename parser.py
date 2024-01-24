
import cmath
import numpy as np

def extract_state_from_file(file, n_qubits):
    with open(file) as f:
        state = []
        while True:
            line = f.readline()
            if len(line) == 0: break
            if len(line) < 3: continue
            pre, amp, phase = line.split()
            amp = float(amp)
            phase = float(phase)
            state.append(amp * cmath.exp(1.0j * phase))
    assert(len(state) == (1 << n_qubits))
    return np.array(state)
