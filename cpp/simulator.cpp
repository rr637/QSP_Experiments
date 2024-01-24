#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <stdexcept>
#include "simulator.hpp"
#include <string>

//int main() {
//    quantum_state state = start_state(2);
//    std::vector<Gate> gates = {
//        Gate::RY(0, 1),
//        Gate::RZ(1, 2),
//        Gate::RX(0, 3),
//        Gate::RX(1, 4),
//        Gate::CX(0, 1),
//        Gate::RY(0, 2),
//        Gate::RY(1, 1),
//        Gate::CX(1, 0),
//    };
//    for(const Gate& g : gates) state = apply_gate(g, state);
//    for (const auto& e : state) {
//        std::cout << e << std::endl;
//    }
//    return 0;
//}

std::vector<double> magnitudes(quantum_state s) {
    std::vector<double> output = std::vector<double>(s.size());
    for(int i = 0; i < s.size(); i++) {
        output[i] = abs(s[i]);
    }
    return output;
}

std::vector<double> phases(quantum_state s) {
    std::vector<double> output = std::vector<double>(s.size());
    for(int i = 0; i < s.size(); i++) {
        output[i] = arg(s[i]);
    }
    return output;
}

double norm(const quantum_state& s) {
    double output = 0;
    for(const std::complex<double>& c : s) output += abs(c) * abs(c);
    return sqrt(output);
}

std::normal_distribution<double> RNG_NORMAL(0, 1);
std::uniform_real_distribution<double> RNG_UNIFORM(-M_PI, M_PI);
auto generator = std::default_random_engine();

quantum_state random_state(int n_qubits) {
    auto amps = std::vector<double>(1 << n_qubits);
    auto phases = std::vector<double>(1 << n_qubits);
    for(int i = 0; i < 1 << n_qubits; i++) {
        amps[i] = RNG_NORMAL(generator);
        phases[i] = RNG_UNIFORM(generator);
    }
    double norm = 0;
    for(const double& a : amps) norm += a * a;
    norm = sqrt(norm);
    quantum_state output = std::vector<std::complex<double> >(1 << n_qubits);
    for(int i = 0; i < 1 << n_qubits; i++) {
        output[i] = std::complex<double>(amps[i] * cos(phases[i]) / norm, 
            amps[i] * sin(phases[i]) / norm);
    }
    return output;
}

quantum_state start_state(int n_qubits) {
    quantum_state output = std::vector<std::complex<double> >(1 << n_qubits, 0);
    output[0] += 1;
    return output;
}

Gate::Gate(GateType type, int target, double angle) : type(type), target(target), angle(angle) {
    
}

Gate Gate::RX(int target, double angle) {
    return Gate(GateType::RX, target, angle);
}

Gate Gate::RY(int target, double angle) {
    return Gate(GateType::RY, target, angle);
}

Gate Gate::RZ(int target, double angle) {
    return Gate(GateType::RZ, target, angle);
}

Gate Gate::CX(int control, int target) {
    if(target == control + 1) {
        return Gate(GateType::CX, control, 0);
    } else if(target == control - 1) {
        return Gate(GateType::XC, target, 0);
    } else {
        throw std::invalid_argument("non-nearest neighbor cx gate");
    }
}

bool Gate::is_cx() const {
    return this->type == GateType::CX || this->type == GateType::XC;
}

Gate Gate::set_angle(double angle) const {
    if(this->is_cx()) throw std::runtime_error("cannot set angle on cx gate");
    return Gate(this->type, this->target, angle);
}

Gate Gate::inverse() const {
    if(this->is_cx()) return *this;
    return this->set_angle(-this->angle);
}

std::string Gate::to_string() const {
    switch(this->type) {
        case GateType::RX: return "rx " + std::to_string(this->target) + " " 
            + std::to_string(this->angle);
        case GateType::RY: return "ry " + std::to_string(this->target) + " "
            + std::to_string(this->angle);
        case GateType::RZ: return "rz " + std::to_string(this->target) + " "
            + std::to_string(this->angle);
        case GateType::CX: return "cx " + std::to_string(this->target) + " "
            + std::to_string(this->target + 1);
        case GateType::XC: return "xc " + std::to_string(this->target + 1) + " "
            + std::to_string(this->target);
        default: throw std::runtime_error("incomplete case match");
    }
}

std::uniform_real_distribution<double> PERTURBATION(-1, 1);
Gate Gate::mutate(double perturbation) const {
    if(this->is_cx()) return *this;
    double delta = perturbation * PERTURBATION(generator);
    Gate output(this->type, this->target, this->angle + delta);
    return output;
}

inline quantum_state apply_rotation_gate(const std::complex<double>& m00, 
    const std::complex<double>& m01, const std::complex<double>& m10,
    const std::complex<double>& m11, const quantum_state& state, int index) {
    quantum_state output = std::vector<std::complex<double> >(state.size());
    for(int i = 0; i < state.size(); i++) {
        if((i >> index) % 2) continue;
        const std::complex<double>& e0 = state[i];
        const std::complex<double>& e1 = state[i ^ (1 << index)];
        output[i] = m00 * e0 + m01 * e1;
        output[i ^ (1 << index)] = m10 * e0 + m11 * e1;
    }
    return output;
}

inline quantum_state apply_cx_gate(int control, int target, const quantum_state& state) {
    quantum_state output = std::vector<std::complex<double> >(state.size());
    for(int i = 0; i < state.size(); i++) {
        if((i >> control) % 2) output[i] = state[i ^ (1 << target)];
        else output[i] = state[i];
    }
    return output;
}

quantum_state apply_gate(const Gate& gate, const quantum_state& state) {
    if(gate.type == GateType::RX) {
        return apply_rotation_gate(
            std::complex<double>(cos(gate.angle / 2), 0),
            std::complex<double>(0, -sin(gate.angle / 2)),
            std::complex<double>(0, -sin(gate.angle / 2)),
            std::complex<double>(cos(gate.angle / 2)),
            state, gate.target);
    }
    if(gate.type == GateType::RY) {
        return apply_rotation_gate(
            std::complex<double>(cos(gate.angle / 2), 0),
            std::complex<double>(-sin(gate.angle / 2), 0),
            std::complex<double>(sin(gate.angle / 2), 0),
            std::complex<double>(cos(gate.angle / 2), 0),
            state, gate.target);
    }
    if(gate.type == GateType::RZ) {
        return apply_rotation_gate(
            std::complex<double>(cos(gate.angle / 2), -sin(gate.angle / 2)),
            0,
            0,
            std::complex<double>(cos(gate.angle / 2), sin(gate.angle / 2)),
            state, gate.target);
    }
    if(gate.type == GateType::CX) {
        return apply_cx_gate(gate.target, gate.target + 1, state);
    }
    if(gate.type == GateType::XC) {
        return apply_cx_gate(gate.target + 1, gate.target, state);
    }
    throw std::runtime_error("incomplete case match");
}

std::string bin_string(int number, int len) {
    std::string output = "";
    for(int i = 0; i < len; i++) {
        if(number % 2) output = "1" + output;
        else output = "0" + output;
        number = number >> 1;
    }
    return output;
}

void print_state(const quantum_state& s, int n_qubits) {
    for(int i = 0; i < s.size(); i++) {
        std::cout 
            << bin_string(i, n_qubits)
            << " "
            << abs(s[i])
            << " "
            << arg(s[i]) - arg(s[0])
            << std::endl;
    }
}

int cx_count(const gate_sequence& gates) {
    int output = 0;
    for(const Gate& g : gates) {
        if(g.is_cx()) output++;
    }
    return output;
}
