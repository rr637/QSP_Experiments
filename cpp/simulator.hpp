
#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <string>

using quantum_state = std::vector<std::complex<double> >;

quantum_state random_state(int n_qubits);
quantum_state start_state(int n_qubits);
std::vector<double> magnitudes(quantum_state);
std::vector<double> phases(quantum_state);
double norm(const quantum_state&);
void print_state(const quantum_state& s, int num_qubits);

enum GateType {
    RX, RY, RZ, CX, XC
};

class Gate {
    Gate(GateType type, int target, double angle);
    public:
        static Gate RX(int target, double angle);
        static Gate RY(int target, double angle);
        static Gate RZ(int target, double angle);
        static Gate CX(int control, int target);
        int target;
        double angle;
        GateType type;
        Gate mutate(double perturbation) const;
        Gate set_angle(double angle) const;
        bool is_cx() const;
        Gate inverse() const;
        std::string to_string() const;
};

quantum_state apply_gate(const Gate&, const quantum_state&);

using gate_sequence = std::vector<Gate>;
int cx_count(const gate_sequence& gates);
