import numpy
import random
import cmath

def rand_state(n):
    state = numpy.array([random.gauss(0, 1) for _ in range(1 << n)])
    norm = state.dot(state) ** 0.5
    return state / norm

def rand_complex_state(n):
    amps = rand_state(n)
    phases = numpy.array([cmath.exp(2.0j * cmath.pi * random.random()) for _ in range(1 << n)])
    return amps * phases

def complex_to_amp_phase(c):
    amp = abs(c)
    phase = cmath.phase(c)
    return amp, phase

def format_number(number):
    # Format to 3 significant figures
    return '{:.6g}'.format(number)

n_qubits = 5
state = rand_complex_state(n_qubits)

# Convert complex numbers to amplitude and phase
formatted_state = [(format(i, f'0{n_qubits}b'), format_number(amp), format_number(phase)) 
                   for i, (amp, phase) in enumerate(map(complex_to_amp_phase, state))]

# Write state to a text file
with open("output.txt", "w") as f:
    for prefix, amp, phase in formatted_state:
        f.write(f"{prefix} {amp} {phase}\n")
