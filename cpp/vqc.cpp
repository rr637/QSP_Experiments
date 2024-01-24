
#include "simulator.hpp"
#include "vqc.hpp"

#include <vector>
#include <complex>
#include <algorithm>
#include <iostream>

using gate_sequence = std::vector<Gate>;

double fidelity(const quantum_state&, const quantum_state&);
//int main() {
//    int n_qubits = 8;
//    int max_cx = 109;
//    for(int j = 0; j < 100; j++) {
//        quantum_state target = random_state(n_qubits);
//        gate_sequence gates;
//        for(int i = 0; i < n_qubits; i++) {
//            gates.push_back(Gate::RY(i, 0));
//            gates.push_back(Gate::RZ(i, 0));
//        }
//        int layer = 0;
//        while(cx_count(gates) < max_cx) {
//            for(int i = layer; i < n_qubits - 1; i += 2) {
//                gates.push_back(Gate::CX(i, i + 1));
//                gates.push_back(Gate::RY(i, 0));
//                gates.push_back(Gate::RZ(i, 0));
//                gates.push_back(Gate::RY(i + 1, 0));
//                gates.push_back(Gate::RZ(i + 1, 0));
//            }
//            layer = (layer + 1) % 2;
//        }
//        while(fitness(target, gates) < 0.95) {
//            gates = gradient_descent(gates, target, 0.1);
//        }
//    }
//}

double fitness(const quantum_state& target, const gate_sequence& seq) {
    quantum_state state = std::vector<std::complex<double> >(target.size());
    state[0] += 1;
    for(const Gate& g : seq) state = apply_gate(g, state);
    return fidelity(state, target);
}

double fidelity(const quantum_state& a, const quantum_state& b) {
    if(a.size() != b.size()) throw std::runtime_error("Unequal dimensions");
    std::complex<double> output = 0;
    for(int i = 0; i < a.size(); i++) {
        output += conj(a[i]) * b[i];
    }
    return abs(output) * abs(output);
}

double compute_gradient(const gate_sequence& gates, const quantum_state& target, int index) {
    const Gate& base = gates[index];
    if(base.is_cx()) return 0;
    quantum_state start = std::vector<std::complex<double> >(target.size());
    start[0] += 1;
    for(int i = 0; i < index; i++) start = apply_gate(gates[i], start);
    quantum_state end = target;
    for(int i = gates.size() - 1; i > index; i--) {
        end = apply_gate(gates[i].inverse(), end);
    }
    //plus case:
    quantum_state plus = apply_gate(base.set_angle(base.angle + M_PI / 2), 
        start);
    double cost_plus = fidelity(plus, end);
    quantum_state minus = apply_gate(base.set_angle(base.angle - M_PI / 2), 
        start);
    double cost_minus = fidelity(minus, end);
    return cost_plus - cost_minus;
}

std::vector<double> compute_gradient(const gate_sequence& gates,
    const quantum_state& start, const quantum_state& target) {
    std::vector<quantum_state> start_cache;
    std::vector<quantum_state> end_cache 
        = std::vector<quantum_state>(gates.size());
    start_cache.reserve(gates.size());
    start_cache.push_back(start);
    for(int i = 0; i < gates.size() - 1; i++) {
        start_cache.push_back(apply_gate(gates[i], start_cache[i]));
    }
    end_cache[gates.size() - 1] = target;
    for(int i = gates.size() - 2; i >= 0; i--) {
        end_cache[i] = apply_gate(gates[i + 1].inverse(), end_cache[i + 1]);
    }
    std::vector<double> output = std::vector<double>(gates.size(), 0);
    for(int i = 0; i < gates.size(); i++) {
        const Gate& base = gates[i];
        if(base.is_cx()) continue;
        quantum_state& init = start_cache[i];
        quantum_state& end = end_cache[i];
        quantum_state plus = apply_gate(base.set_angle(base.angle + M_PI / 2),
            init);
        double cost_plus = fidelity(plus, end);
        quantum_state minus = apply_gate(base.set_angle(base.angle - M_PI / 2),
            init);
        double cost_minus = fidelity(minus, end);
        output[i] = cost_plus - cost_minus;
    }
    return output;
}

std::vector<double> compute_gradient(const gate_sequence& gates, 
    const quantum_state& target) {
    //start_cache[i] = result of applying the first [i] gates in [gates] to
    //    state |0>, not including the gate at index [i].
    //end_cache[i] = result of applying the inverse of every gate after index 
    //    [i] in [gates], in reverse order, not including the gate at [i], to
    //    the target state.
    quantum_state start = std::vector<std::complex<double> >(target.size(), 0);
    start[0] += 1;
    return compute_gradient(gates, start, target);
}

gate_sequence gradient_descent(const gate_sequence& gates, const quantum_state& target,
    double step_size) {
    gate_sequence output;
    output.reserve(gates.size());
    std::vector<double> gradient = compute_gradient(gates, target);
    for(int i = 0; i < gates.size(); i++) {
        if(gradient[i] == 0) output.push_back(gates[i]);
        else {
            output.push_back(gates[i].set_angle(gates[i].angle + 
                step_size * gradient[i]));
        }
    }
    return output;
}

gate_sequence gradient_descent(const gate_sequence& gates, const quantum_state& start, const quantum_state& target, double step_size) {
    gate_sequence output;
    output.reserve(gates.size());
    std::vector<double> gradient = compute_gradient(gates, start, target);
    for(int i = 0; i < gates.size(); i++) {
        if(gradient[i] == 0) output.push_back(gates[i]);
        else {
            output.push_back(gates[i].set_angle(gates[i].angle + step_size * gradient[i]));
        }
    }
    return output;
}

