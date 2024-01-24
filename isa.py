
from qc import Gate, apply_gate, substate, rand_state, print_state
import cmath
from math import atan, acos
from util import bit_flip, max_index, norm
import numpy as np
from sp3_extension import sp3_1, sp2

def templates(n_qubits):
    assert(n_qubits > 3)
    output = []
    for i in range(1, n_qubits - 1):
        output.append("?" * i + "**" + "0" * (n_qubits - i - 2))
        output.append("0" * (n_qubits - i - 2) + "**" + "?" * i)
    return output

def str_replace(string, char, index):
    assert(index >= 0 and index < len(string))
    output = string[0:index]
    output += char
    output += string[index + 1:]
    return output

def template_to_pattern(template):
    index = template.find("?")
    if index == -1: return [template]
    zero_branch = str_replace(template, "0", index)
    one_branch = str_replace(template, "1", index)
    output = template_to_pattern(zero_branch)
    output.extend(template_to_pattern(one_branch))
    return output

def patterns(n_qubits):
    t = templates(n_qubits)
    output = []
    for tt in t: output.extend(template_to_pattern(tt))
    for i in range(n_qubits):
        output.append("0" * i + "*" + "0" * (n_qubits - i - 1))
    return output

#The number of nnCX's required to transform branch to a string with all zeroes,
# except one 1 on the far left.
def branch_cx_cost(branch):
    i1 = branch.rfind("1")
    if i1 == -1: return -2
    #output is number of zeroes between i1 and left plus number of chars to the 
    # left of i1
    return sum(int(branch[i] == "0") for i in range(i1)) + i1

def pattern_cx_cost(pattern):
    star_index = pattern.index("*")
    if star_index == len(pattern) - 1 or pattern[star_index + 1] != "*": 
        return 0
    left_branch = pattern[0:star_index]
    right_branch = pattern[star_index + 2:]
    if left_branch == "0" * len(left_branch):
        return 3 + branch_cx_cost(right_branch)
    return 3 + branch_cx_cost(left_branch[::-1])

def find_target(src, dest):
    xor = int(src) ^ int(dest)
    output = xor.bit_length() - 1
    assert(1 << output == xor)
    return output

def find_control_target(src, dest):
    target = find_target(src, dest)
    if (src >> (target + 1)) % 2: return target + 1, target
    return target - 1, target

class StateTracker:
    def __init__(self, state):
        self.state = state
        self.gates = [[]] #Invariant: self.gates always has at least one block
        self.n_qubits = len(state).bit_length() - 1
    def new_block(self):
        self.gates.append([])
    def apply_gate(self, gate):
        self.state = apply_gate(gate, self.state)
        self.gates[-1].append(gate)
    def apply_ry(self, target, angle):
        self.apply_gate(Gate.RY(target, angle, self.n_qubits))
    def apply_rz(self, target, angle):
        self.apply_gate(Gate.RZ(target, angle, self.n_qubits))
    def apply_cx(self, control, target):
        self.apply_gate(Gate.CX(control, target, self.n_qubits))
    #If there are no gates in the current block, throws an error
    def undo_gate(self):
        if len(self.gates[-1]) == 0: raise
        gate = self.gates[-1][-1]
        self.gates[-1].pop()
        self.state = apply_gate(gate.inverse(), self.state)
    def undo_block(self):
        block = self.gates[-1][:]
        if len(self.gates) == 1: self.gates[0].clear()
        else: self.gates.pop()
        for gate in reversed(block):
            self.state = apply_gate(gate.inverse(), self.state)
    def rotate_merge(self, src, dest):
        self.unify_phase(src, dest)
        target = find_target(src, dest)
        angle = -2 * atan(abs(self.state[src] / self.state[dest]))
        if dest > src: angle *= -1
        self.apply_ry(target, angle)
    def control_rotate_merge(self, src, dest):
        self.unify_phase(src, dest)
        control, target = find_control_target(src, dest)
        angle = atan(abs(self.state[dest] / self.state[src]))
        if dest > src: angle *= -1
        self.apply_ry(target, angle)
        self.apply_cx(control, target)
        self.apply_ry(target, -angle)
    def unify_phase(self, src, dest):
        target = find_target(src, dest)
        dphase = cmath.phase(self.state[dest]) - cmath.phase(self.state[src])
        if dest > src: dphase *= -1
        self.apply_rz(target, dphase)
    def _rotate_chunk_angles(self, src_pattern, target):
        dest_pattern = None
        down_rotate = False #true if dest < src
        v0 = None
        v1 = None
        if src_pattern[target] == "0": 
            dest_pattern = str_replace(src_pattern, "1", target)
            v0 = substate(self.state, src_pattern)
            v1 = substate(self.state, dest_pattern)
        elif src_pattern[target] == "1":
            dest_pattern = str_replace(src_pattern, "0", target)
            v0 = substate(self.state, dest_pattern)
            v1 = substate(self.state, src_pattern)
            down_rotate = True
        else:
            raise
        v0_mag = sum(abs(v0) * abs(v0))
        v1_mag = sum(abs(v1) * abs(v1))
        cross = sum(vv0.conjugate() * vv1 for vv0, vv1 in zip(v0, v1))
        cross_mag = abs(cross)
        phi = -cmath.phase(cross)
        A = (v0_mag - v1_mag) / 2
        theta = -acos(A / (A * A + cross_mag * cross_mag) ** 0.5)
        if not down_rotate: theta += cmath.pi
        return phi, theta
    #target refers to the index in the pattern
    def rotate_merge_chunk(self, src_pattern, target):
        phi, theta = self._rotate_chunk_angles(src_pattern, target)
        target = self.n_qubits - 1 - target
        self.apply_rz(target, phi)
        self.apply_ry(target, theta)
    #control and target refer to indices in the pattern
    def control_rotate_merge_chunk(self, src_pattern, control, target):
        phi, theta = self._rotate_chunk_angles(src_pattern, target)
        theta = (theta - cmath.pi) / 2
        target = self.n_qubits - 1 - target
        control = self.n_qubits - 1 - control
        self.apply_rz(target, phi)
        self.apply_ry(target, theta)
        self.apply_cx(control, target)
        self.apply_ry(target, -theta)
    def gate_sequence(self):
        output = []
        for block in self.gates:
            output.extend(block)
        return output

def pattern_value(pattern, state):
    if pattern_cx_cost(pattern) <= 2: 
        return norm(substate(state, pattern)) ** 2 - abs(state[0]) ** 2
    return norm(substate(state, pattern)) ** 2

def list_cx_targets(pattern):
    #list everything next to a 1, as long as it's not a *
    output = set()
    for i in range(len(pattern) - 1):
        if pattern[i] == "1" and pattern[i + 1] != "*": output.add(i + 1)
        if pattern[i + 1] == "1" and pattern[i] != "*": output.add(i)
    return list(output)

def pattern_bit_flip(pattern, target):
    if pattern[target] == "0": return str_replace(pattern, "1", target)
    if pattern[target] == "1": return str_replace(pattern, "0", target)
    raise

def pattern_target_value(pattern, target, state):
    subs1 = substate(state, pattern)
    subs2 = substate(state, pattern_bit_flip(pattern, target))
    subs1_mag = np.vdot(subs1, subs1)
    subs2_mag = np.vdot(subs2, subs2)
    subs12_mag = abs(np.vdot(subs1, subs2))
    A = 0.5 * (subs1_mag - subs2_mag)
    output = 0.5 * (subs1_mag + subs2_mag) + \
        (A * A + subs12_mag * subs12_mag) ** 0.5
    return output

def pattern_merge_cost(pattern, target):
    cx_cost = pattern_cx_cost(pattern)
    if pattern_cx_cost(pattern_bit_flip(pattern, target)) > cx_cost:
        return 1 + cx_cost
    return cx_cost

def isa_prepare(state, target_fidelity=0.95):
    assert(abs(1 - np.vdot(state, state)) < 0.001)
    tracker = StateTracker(state)
    #Phase 1: same as original greedy
    phase_1_candy = [i for i in range(tracker.n_qubits)]
    index = max_index(abs(tracker.state))
    while len(phase_1_candy) > 0:
        values = [abs(tracker.state[bit_flip(index, m)]) for m in phase_1_candy]
        target = phase_1_candy[max_index(values)]
        src = max(index, bit_flip(index, target))
        dest = min(index, bit_flip(index, target))
        tracker.rotate_merge(src, dest)
        phase_1_candy.remove(target)
        index = dest

    #Phase 2: iteratively select a pattern then implement that pattern
    pats = patterns(tracker.n_qubits)
    costs = np.array([pattern_cx_cost(p) for p in pats])
    past_pattern = None
    while abs(tracker.state[0]) ** 2 < target_fidelity:
        values = np.array([pattern_value(p, tracker.state) for p in pats])
        value_per_cost = values / (costs + 1)
        pattern = pats[max_index(value_per_cost)]
        if past_pattern is not None \
          and past_pattern.find("*") == pattern.find("*"):
            tracker.undo_block()
        while pattern_cx_cost(pattern) > 3:
            #Make list of indices than can be flipped
            targets = list_cx_targets(pattern)
            values = np.array([pattern_target_value(pattern, target, \
              tracker.state) for target in targets])
            merge_costs = np.array([pattern_merge_cost(pattern, target) \
              for target in targets])
            target = targets[max_index(values / merge_costs)]
            if target == 0: control = 1
            elif pattern[target - 1] == "1": control = target - 1
            else: control = target + 1
            if pattern_cx_cost(pattern_bit_flip(pattern, target)) < \
              pattern_cx_cost(pattern):
                src = pattern
                pattern = pattern_bit_flip(pattern, target)
            else:
                src = pattern_bit_flip(pattern, target)
            tracker.control_rotate_merge_chunk(src, control, target)
        #now use the sp3 extension to merge.
        #find stars, figure out which side is the 1. Then, compute q0, q1, q2
        star_index = pattern.find("*")
        if pattern_cx_cost(pattern) == 3:
            if star_index == -1: raise
            if star_index == 0:
                q0 = tracker.n_qubits - 1
                q1 = tracker.n_qubits - 2
                q2 = tracker.n_qubits - 3
            elif pattern[star_index - 1] == "1":
                q0 = tracker.n_qubits - 2 - star_index
                q1 = tracker.n_qubits - 1 - star_index
                q2 = tracker.n_qubits - star_index
            elif pattern[star_index + 2] == "1":
                q0 = tracker.n_qubits - 1 - star_index
                q1 = tracker.n_qubits - 2 - star_index
                q2 = tracker.n_qubits - 3 - star_index
            tracker.new_block()
            sp3_1(tracker, q0, q1, q2)
            tracker.new_block()
            sp2(tracker, q0, q1)
        elif pattern_cx_cost(pattern) == 1:
            q0 = tracker.n_qubits - 1 - star_index
            q1 = tracker.n_qubits - 2 - star_index
            sp2(tracker, q0, q1)
            past_pattern = None
        else:
            target = tracker.n_qubits - 1 - star_index
            tracker.rotate_merge((1 << target), 0)
            past_pattern = None
    return tracker.gate_sequence()

