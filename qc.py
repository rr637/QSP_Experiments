
import math
import cmath
import qiskit
import numpy as np
import random
from util import get_bit

random.seed(42)

#Class representing quantum gates. Use the static methods to construct new
# objects. 
class Gate:
    #cx is a CX gate with the control on the lower index qubit; xc is a CX
    # gate with the control on the higher index qubit.
    types = ["rx", "ry", "rz", "cx", "xc"]
    def __init__(self, gate_type, target, angle, n_qubits):
        self.gate_type = gate_type
        self.target = target
        self.angle = angle
        self.n_qubits = n_qubits

    def is_cx(self):
        return self.gate_type in ["cx", "xc"]
    
    def inverse(self):
        if self.is_cx(): return self
        return Gate(self.gate_type, self.target, -self.angle, self.n_qubits)

    def to_string(self):
        if self.is_cx():
            control = self.target
            target = control + 1
            if self.gate_type == "xc": 
                control -= 1
                target += 1
            return "CX control=" + str(control) + " target=" + str(target)
        return self.gate_type.upper() + " target=" + str(self.target) + \
        " angle=" + str(self.angle)

    #CX requires control and target to be on nearest neighbor qubits
    @staticmethod
    def CX(control, target, n_qubits):
        if target == control + 1:
            return Gate("cx", control, 0, n_qubits)
        if control == target + 1:
            return Gate("xc", target, 0, n_qubits)
        raise RuntimeError()

    @staticmethod
    def RY(target, angle, n_qubits):
        return Gate("ry", target, angle, n_qubits)

    @staticmethod
    def RX(target, angle, n_qubits):
        return Gate("rx", target, angle, n_qubits)

    @staticmethod
    def RZ(target, angle, n_qubits):
        return Gate("rz", target, angle, n_qubits)

#Helper function
def apply_rotation_gate(gate, state, index):
    output = [None for _ in range(len(state))]
    for i in range(len(output)):
        if output[i] is not None: continue
        e0 = state[i]
        # print(f"Type of index: {type(index)}")
        # print(f"Value of index: {index}")
        e1 = state[i + (1 << index)]

        output[i] = gate[0][0] * e0 + gate[0][1] * e1
        output[i + (1 << index)] = gate[1][0] * e0 + gate[1][1] * e1
    return np.array(output)

#Helper function
def apply_cx_gate(state, control, target):
    output = [None for _ in range(len(state))]
    for i in range(len(output)):
        if (i >> control) % 2: output[i] = state[i ^ (1 << target)]
        else: output[i] = state[i]
    return np.array(output)

#Given Gate object [gate] and list of complex amplitudes [state] representing a
# quantum state, returns a new list of complex amplitudes representing the
# result of applying quantum gate [gate] to quantum state [state]. Does not
# modify [state].
def apply_gate(gate, state):
    if gate.gate_type == "rx":
        mat = [[math.cos(gate.angle / 2), -1.0j * math.sin(gate.angle / 2)], \
               [-1.0j * math.sin(gate.angle / 2), math.cos(gate.angle / 2)]]
        return apply_rotation_gate(mat, state, gate.target)
    if gate.gate_type == "ry":
        mat = [[math.cos(gate.angle / 2), -math.sin(gate.angle / 2)], \
               [math.sin(gate.angle / 2), math.cos(gate.angle / 2)]]
        return apply_rotation_gate(mat, state, gate.target)
    if gate.gate_type == "rz":
        mat = [[cmath.exp(-1.0j * gate.angle / 2), 0], \
               [0, cmath.exp(1.0j * gate.angle / 2)]]
        return apply_rotation_gate(mat, state, gate.target)
    if gate.gate_type == "cx":
        return apply_cx_gate(state, gate.target, gate.target + 1)
    if gate.gate_type == "xc":
        return apply_cx_gate(state, gate.target + 1, gate.target)

#Helper function
def validate_apply_gate():
    circuit = qiskit.QuantumCircuit(2)
    circuit.ry(1, 0)
    circuit.rz(2, 1)
    circuit.rx(3, 0)
    circuit.rx(4, 1)
    circuit.cx(0, 1)
    circuit.ry(2, 0)
    circuit.ry(1, 1)
    circuit.cx(1, 0)
    simulator = qiskit.Aer.get_backend("statevector_simulator")
    sv1 = np.array(simulator.run(circuit).result().get_statevector())
    sv2 = [1, 0, 0, 0]
    for gate in [Gate.RY(0, 1, 2), Gate.RZ(1, 2, 2), Gate.RX(0, 3, 2), 
    Gate.RX(1, 4, 2), Gate.CX(0, 1, 2), Gate.RY(0, 2, 2), Gate.RY(1, 1, 2), 
    Gate.CX(1, 0, 2)]:
        sv2 = apply_gate(gate, sv2)
    sv2 = np.array(sv2)
    delta = sv1 - sv2
    assert(delta.dot(delta) < 0.00001)
    print("Success")

#Returns a random complex-amplitude quantum state on [n] qubits
def rand_state(n):
    amps = np.array([random.gauss(0, 1) for _ in range(1 << n)])
    norm = sum(amps * amps)
    amps /= (norm ** 0.5)
    phases = np.array([cmath.exp(cmath.pi * 2.0j * random.random()) \
        for _ in range(1 << n)])
    return amps * phases

def bin_string(n, bits=None):
    if bits is None: bits = n.bit_length()
    buf = []
    for i in range(bits - 1, -1, -1):
        buf.append(str(get_bit(n, i)))
    return "".join(buf)

def print_state(state):
    bits = len(state).bit_length() - 1
    for i, a in enumerate(state):
        if abs(a) < 1e-7:
            amp = " 0.0000000"
            phase = "  0.000000"
        else:
            amp = " " + str(abs(a))[0:9]
            phase = cmath.phase(a) - cmath.phase(state[0])
            if abs(phase) < 1e-6: phase = "  0.000000"
            elif phase < 0: phase = " " + str(phase)[0:9]
            else: phase = "  " + str(phase)[0:8]
        print(bin_string(i, bits), amp, phase)

def matches(string, pattern):
    if len(string) != len(pattern): return False
    for p, c in zip(pattern, string):
        if p != "*" and p != c: return False
    return True

def substate(state, pattern):
    if type(pattern) is not str: raise
    bits = len(pattern)
    output = []
    for i, a in enumerate(state):
        string = bin_string(i, bits)
        if matches(string, pattern): output.append(a)
    return np.array(output)

if __name__ == "__main__":
    validate_apply_gate()

