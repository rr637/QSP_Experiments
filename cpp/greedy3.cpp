
#include "greedy3.hpp"
#include "sp3.hpp"
#include "state_tracker.hpp"
#include "pattern.hpp"
#include "simulator.hpp"
#include "greedy.hpp"
#include <vector>
#include <iostream>
#include <algorithm>

gate_sequence greedy3_prepare(const quantum_state& state, int n_qubits,
  double target_fidelity);

//int main() {
//    int n_qubits = 8;
//    int count = 100;
//    double total = 0.0;
//    for(int i = 0; i < count; i++) {
//        quantum_state state = random_state(n_qubits);
//        gate_sequence gates = greedy3_prepare(state, n_qubits, 0.95);
//        for(const Gate& g : gates) state = apply_gate(g, state);
//        if(abs(state[0]) * abs(state[0]) < 0.94999) {
//            throw std::runtime_error("Infidelity");
//        }
//        for(const Gate& g : gates) {
//            if(g.is_cx()) total += 1;
//        }
//    }
//    std::cout << total / 100 << std::endl;
//}

std::vector<Pattern> generate_patterns(int n_qubits) {
    std::vector<Pattern> output;
    //first generate all the lower patterns ...
    //00 ... 0**
    //00 ... 0**#
    //...
    //**##...##
    for(int k = 1; k < n_qubits - 1; k++) {
        for(int j = 0; j < (1 << k); j++) {
            output.emplace_back(j, 0, k, k + 2);
        }
    }
    //now generate all the upper patterns ...
    //**0 ... 0
    //#**0... 0
    //...
    //## ... #**
    for(int k = 1; k < n_qubits - 1; k++) {
        for(int j = 1; j < (1 << k); j++) {
            output.emplace_back(0, j, n_qubits - 2 - k, n_qubits - k);
        }
    }
    output.emplace_back(0, 0, 0, 2);
    //now generate all the patterns with upper = lower = 0, but one star in 
    //the middle
    for(int k = 0; k < n_qubits; k++) {
        output.emplace_back(0, 0, k, k + 1);
    }
    return output;
}

int pattern_cost(const Pattern& p) {
    if(p.middle_end - p.middle_start == 1) return 0;
    if(p.middle_end - p.middle_start != 2) {
        throw std::runtime_error("Unimplemented");
    }
    if(p.lower == 0 && p.upper == 0) return 1;
    if(p.lower == 0) {
        int u = p.upper;
        int output = 3;
        while(u > 1) {
            if(u & 1) output += 1;
            else output += 2;
            u = u >> 1;
        }
        return output;
    }
    if(p.upper == 0) {
        int l = p.lower;
        int output = 2;
        int bits = p.middle_start;
        while(!(l & 1)) {
            l = l >> 1;
            bits--;
        }
        for(int i = 0; i < bits; i++) {
            if(l & 1) output += 1;
            else output += 2;
            l = l >> 1;
        }
        return output;
    }
    throw std::runtime_error("Lower and upper are both nonzero.");
}

std::vector<int> list_cx_targets(const Pattern& p, int n_qubits) {
    std::vector<int> output;
    if(p.lower == 0) {
        for(int i = p.middle_end; i < n_qubits - 1; i++) {
            if((p.upper & (1 << (i - p.middle_end + 1))) 
              && (output.size() == 0 || output.back() != i)) {
                output.push_back(i);
            }
            if(p.upper & (1 << (i - p.middle_end))) output.push_back(i + 1);
        }
        return output;
    }
    //then p.upper is 0
    for(int i = 0; i < p.middle_start - 1; i++) {
        if((p.lower & (1 << (i + 1))) 
          && (output.size() == 0 || output.back() != i)) {
            output.push_back(i);
        }
        if(p.lower & (1 << i)) output.push_back(i + 1);
    }
    return output;
}

gate_sequence greedy3_prepare(const quantum_state& state, int n_qubits, 
  double target_fidelity) {
    StateTracker tracker(state);
    greedy_phase1(tracker, n_qubits);
    std::vector<Pattern> patterns = generate_patterns(n_qubits);
    std::vector<int> costs;
    for(const Pattern& p : patterns) costs.push_back(pattern_cost(p));
    int last2 = -1;
    while(abs(tracker.state[0]) * abs(tracker.state[0]) < target_fidelity) {
        std::vector<double> values;
        for(int i = 0; i < patterns.size(); i++) {
            Pattern& p = patterns[i];
            double discount = 0;
            if(p.lower == p.upper) {
                discount = abs(tracker.state[0]) * abs(tracker.state[0]);
            }
            double subs_norm = norm(substate(tracker.state, p));
            double increase = subs_norm * subs_norm - discount;
            values.push_back(increase / (costs[i] + 1));
        }
        Pattern p_selected = patterns[max_index(values)];
        if(p_selected.middle_end - p_selected.middle_start == 2 
          && p_selected.middle_start == last2) {
            int old = tracker.cx_count();
            tracker.undo_block();
            int new_count = tracker.cx_count();
        }
        tracker.new_block();
        while(pattern_cost(p_selected) > 3) {
            std::vector<int> targets = list_cx_targets(p_selected, n_qubits);
            std::vector<double> t_values;
            quantum_state subs1 = substate(tracker.state, p_selected);
            double subs1_norm = norm(subs1);
            double subs1_mag = subs1_norm * subs1_norm;
            for(int t : targets) {
                Pattern other = p_selected.bit_flip(t);
                quantum_state subs2 = substate(tracker.state, other);
                double subs2_norm = norm(subs2);
                double subs2_mag = subs2_norm * subs2_norm;
                std::complex<double> cross = 0;
                for(int i = 0; i < subs1.size(); i++) {
                    cross += conj(subs1[i]) * subs2[i];
                }
                double subs12_mag = abs(cross);
                double A = 0.5 * (subs1_mag - subs2_mag);
                double value = 0.5 * (subs1_mag + subs2_mag) 
                    + sqrt(A * A + subs12_mag * subs12_mag);
                int cost = std::max(pattern_cost(other), 
                  pattern_cost(p_selected));
                t_values.push_back(value / cost);
            }
            int target = targets[max_index(t_values)];
            Pattern src = p_selected;
            Pattern dest = p_selected.bit_flip(target);
            if(pattern_cost(src) < pattern_cost(dest)) {
                src = p_selected.bit_flip(target);
                dest = p_selected;
            }
            p_selected = dest;
            tracker.control_rotate_merge_chunk(src, dest);
        }
        if(pattern_cost(p_selected) == 3) {
            int q0, q1, q2;
            if(p_selected.lower == 0) {
                q0 = p_selected.middle_start;
                q1 = q0 + 1;
                q2 = q0 + 2;
                last2 = q0;
            } else {
                q0 = p_selected.middle_start + 1;
                q1 = q0 - 1;
                q2 = q0 - 2;
                last2 = q1;
            }
            tracker.new_block();
            sp3_1(tracker, q0, q1, q2);
            tracker.new_block();
            sp2(tracker, q0, q1);
        } else if(pattern_cost(p_selected) == 1) {
            int q0 = p_selected.middle_start;
            sp2(tracker, q0 + 1, q0);
            last2 = q0;
        } else {
            int target = p_selected.middle_start;
            tracker.rotate_merge((1 << target), 0);
            last2 = -1;
        }
    }
    return tracker.gates_seq();
}

