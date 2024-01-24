import numpy as np
import cmath
from math import atan, sin, cos
from qc import Gate, apply_gate, print_state, rand_state
from util import get_bit

def normalize(v):
    norm = sum(abs(v) * abs(v))
    return v / (norm ** 0.5)

def phase_factor(p):
    return cmath.exp(1.0j * p)

def eig_x(abt):
    v1 = normalize(np.array([-abt[0][1], abt[0][0] - 1]))
    v2 = normalize(np.array([-abt[0][1], abt[0][0] + 1]))
    p = np.array([[v1[0], v2[0]], [v1[1], v2[1]]])
    u = p / cmath.sqrt(2) @ np.array([[1, 1], [1, -1]])
    return u

#applies the single-qubit rotation [mat] to the [index]th qubit on tracker's
# state
def zyz(tracker, mat, index):
    global_phase = cmath.phase(mat[0][0])
    z1 = cmath.phase(mat[1][0]) - global_phase
    y = 2 * atan(abs(mat[1][0]) / abs(mat[0][0]))
    z2 = cmath.phase(mat[0][1]) - global_phase + cmath.pi
    tracker.apply_rz(index, z2)
    tracker.apply_ry(index, y)
    tracker.apply_rz(index, z1)

def substate(state, *args):
    if len(args) < 1: raise
    output = []
    for i in range(1 << len(args)):
        index = sum((1 << n) * ((i >> j) & 1) for j, n in enumerate(args))
        output.append(state[index])
    return output

class SubstateView:
    def __init__(self, tracker, *args):
        if len(args) < 1: raise
        self.tracker = tracker
        self.indices = args
        self.n_qubits = len(args)
    def __len__(self):
        return 1 << self.n_qubits
    def __getitem__(self, n):
        if n < 0 or n >= len(self): raise
        index = sum((1 << k) * ((n >> j) & 1) for j, k in enumerate(self.indices))
        return self.tracker.state[index]

#applies a sequence of gates so that the substate of q0, q1 has q1 = 0
def sp2_1(tracker, q0, q1):
    tracker.rotate_merge((1 << q1), 0)
    tracker.control_rotate_merge((1 << q0) + (1 << q1), (1 << q0))
    #tracker.rotate_merge((1 << q0), 0)

#Returns a / b. If b is zero and default is specified, then returns default
# otherwise returns (a + 1e-9) / (b + 1e-9)
def safe_div(a, b, default=None):
    if b != 0: return a / b
    if default is not None: return default
    return (a + 1e-9) / (b + 1e-9)

def safe_option_div(a, b, c, d):
    if abs(a) + abs(b) > abs(c) + abs(d): return safe_div(a, b)
    return safe_div(c, d)

#applies a sequence of gates to zero out q2.
def sp3_1(tracker, q0, q1, q2):
    state = SubstateView(tracker, q0, q1, q2)
    #step d1: setting up the ratios
    if abs(state[7] * state[1] - state[5] * state[3]) < 1e-9:
        r = -safe_option_div(state[3], state[1], state[7], state[5]).conjugate()
        c = safe_option_div(state[5], state[1], state[7], state[3])
        l = safe_div(state[4] - c * state[0], state[6] - c * state[2])
    elif abs(state[6] * state[0] - state[4] * state[2]) < 1e-9:
        l = -safe_option_div(state[2], state[0], state[6], state[4]).conjugate()
        c = safe_option_div(state[4], state[0], state[6], state[2])
        r = safe_div(state[5] - c * state[1], state[7] - c * state[3])
    else:
        A = state[4] * state[1] - state[0] * state[5]
        B = state[4] * state[3] - state[0] * state[7]
        C = state[6] * state[1] - state[2] * state[5]
        D = state[6] * state[3] - state[2] * state[7]
        A_c = A.conjugate()
        B_c = B.conjugate()
        C_c = C.conjugate()
        D_c = D.conjugate()
        
        a = -(C_c * D + A_c * B)
        b = A_c * A - B_c * B + C_c * C - D_c * D
        c = C * D_c + A * B_c
    
        r = safe_div(-b + cmath.sqrt(b * b - 4 * a * c), 2 * a)
        l = safe_div(-(C_c * r + D_c), A_c * r + B_c)
    tr = atan(abs(r))
    pr = safe_div(r, abs(r))
    tl = atan(abs(l))
    pl = safe_div(l, abs(l))

    ml = np.array([[cos(tl), -sin(tl) * pl], [sin(tl), cos(tl) * pl]])
    mr = np.array([[cos(tr), -sin(tr) * pr], [sin(tr), cos(tr) * pr]])

    abt = ml @ mr.transpose().conj()
    global_phase = cmath.phase(abt[0][0])
    gamma = cmath.pi - cmath.phase(abt[1][1]) + global_phase

    ml = np.diag([1, phase_factor(gamma)]) @ ml / phase_factor(global_phase)
    abt = ml @ mr.transpose().conj()
    u = eig_x(abt)
    v = u.transpose().conj() @ ml
    zyz(tracker, v, q1)
    tracker.apply_cx(q0, q1)
    zyz(tracker, u, q1)
    
    #sp2 to zero out q2, but accounting for possible 0/0
    if abs(state[4]) + abs(state[0]) > abs(state[5]) + abs(state[1]):
        tracker.rotate_merge((1 << q2), 0)
    else: 
        tracker.rotate_merge((1 << q2) + (1 << q0), (1 << q0))
    if abs(state[6]) + abs(state[2]) > abs(state[7]) + abs(state[3]):
        tracker.control_rotate_merge((1 << q2) + (1 << q1), (1 << q1))
    else:
        tracker.control_rotate_merge((1 << q2) + (1 << q1) + (1 << q0), \
            (1 << q1) + (1 << q0))
    if abs(state[2]) + abs(state[0]) > abs(state[3]) + abs(state[1]):
        tracker.rotate_merge((1 << q1), 0)
    else:
        tracker.rotate_merge((1 << q1) + (1 << q0), (1 << q0))

def sp2(tracker, q0, q1):
    sp2_1(tracker, q0, q1)
    tracker.rotate_merge((1 << q0), 0)

def sp3(tracker, q0, q1, q2):
    sp3_1(tracker, q0, q1, q2) #zero out q2
    sp2_1(tracker, q0, q1) #zero out q1
    tracker.rotate_merge((1 << q0), 0) #step sp2-2: zero out q0

if __name__ == "__main__":
#    from greedy3 import StateTracker
#    masks = []
#    for i in range(1, 64):
#        masks.append(np.array([get_bit(i, j) for j in range(8)]))
#        print(masks[-1])
#    for mask in masks:
#        state = rand_state(3)
#        state *= mask
#        state = normalize(state)
#        tracker = StateTracker(state)
#        sp3(tracker, 0, 1, 2)
#        assert(abs(1 - abs(tracker.state[0])) >= 1e-6)
#    print("Success")
    from greedy3 import StateTracker
    amps = [
        0.371019,
        0.121945,
        0.413835,
        0.0156576,
        0.270587,
        0.559836,
        0.496921,
        0.205924,
    ]
    phases = [
        -1.61156,
        -2.76215,
        3.11071,
        2.93181,
        -0.241695,
        2.2493,
        -3.04802,
        -1.32082,
    ]
    state = np.array([a * cmath.exp(p * 1.0j) for a, p in zip(amps, phases)])
    print_state(state)
    tracker = StateTracker(state)
    sp3(tracker, 0, 1, 2)
