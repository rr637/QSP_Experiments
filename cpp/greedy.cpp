
#include "greedy.hpp"
#include "simulator.hpp"
#include "state_tracker.hpp"
#include <algorithm>
#include <vector>
#include <iostream>
#include <complex>

//int main() {
//    int n_qubits = 14;
//    int count = 100;
//    double total = 0.0;
//    for(int i = 0; i < count; i++) {
//        std::cout << i << std::endl;
//        quantum_state state = random_state(n_qubits);
//        gate_sequence gates = prepare_greedy(state, n_qubits);
//        total += cx_count(gates);
//    }
//    std::cout << total / count << std::endl;
//}

int top_one(int number) {
    int output = 0;
    while(number > 1) {
        output += 1;
        number = number >> 1;
    }
    return output;
}

int bot_one(int number) {
    int output = 0;
    while(number % 2 == 0) {
        output += 1;
        number = number >> 1;
    }
    return output;
}

bool get_bit(int number, int position) {
    return (number >> position) % 2;
}

int max_index(std::vector<double> list) {
    int output = 0;
    double value = list[0];
    for(int i = 1; i < list.size(); i++) {
        if(list[i] > value) {
            value = list[i];
            output = i;
        }
    }
    return output;
}

void print_list(std::vector<double> list) {
    std::cout << "[";
    for(double d : list) std::cout << d << ", ";
    std::cout << "]" << std::endl;
}

void print_list(std::vector<int> list) {
    std::cout << "[";
    for(int d : list) std::cout << d << ", ";
    std::cout << "]" << std::endl;
}

int cx_dist(int number) {
    int top = top_one(number);
    int bot = bot_one(number);
    int zeroes = 0;
    for(int i = bot + 1; i < top; i++) {
        if((number >> i) % 2) continue;
        zeroes++;
    }
    return zeroes + top - bot;
}

int apply_cx(int number, int control, int target) {
    if((number >> control) % 2) return number ^ (1 << target);
    return number;
}

std::vector<std::pair<int, int> > cflips(int number, int len) {
    std::vector<std::pair<int, int> > output;
    for(int i = 0; i < len - 1; i++) {
        if((number >> i) % 2) output.push_back(std::pair<int, int>(i, i + 1));
        if((number >> (i + 1)) % 2) {
            output.push_back(std::pair<int, int>(i + 1, i));
        }
    }
    return output;
}

void greedy_phase1(StateTracker& tracker, int n_qubits) {
    std::vector<int> candy1;
    for(int i = 0; i < n_qubits; i++) candy1.push_back(i);
    std::vector<double> amps = magnitudes(tracker.state);
    int index = max_index(amps);
    while(candy1.size() > 0) {
        std::vector<double> values;
        for(int c : candy1) {
            values.push_back(abs(tracker.state[index ^ (1 << c)]));
        }
        int target = candy1[max_index(values)];
        int src = std::max(index, index ^ (1 << target));
        int dest = std::min(index, index ^ (1 << target));
        tracker.rotate_merge(src, dest);
        candy1.erase(std::remove(candy1.begin(), candy1.end(), target), 
            candy1.end());
        index = dest;
    }
}

void greedy_phase2_iter(StateTracker& tracker, int n_qubits) {
    std::vector<double> values;
    for(int i = 0; i < 1 << n_qubits; i++) {
        if(i == 0) {
            values.push_back(0);
            continue;
        }
        values.push_back(abs(tracker.state[i]) * abs(tracker.state[i]) / (1 + cx_dist(i)));
    }
    int index = max_index(values);
    while(cx_dist(index) > 0) {
        std::vector<std::pair<int, int> > candy = cflips(index, n_qubits);
        std::vector<double> values;
        for(std::pair<int, int>& cx : candy) {
            double v = abs(tracker.state[index]) * abs(tracker.state[index]);
            int other = index ^ (1 << cx.second);
            v += abs(tracker.state[other]) * abs(tracker.state[other]);
            double denom = 1 + cx_dist(index);
            if(cx_dist(other) > cx_dist(index)) denom += 1;
            values.push_back(v / denom);
        }
        std::pair<int, int> cx = candy[max_index(values)];
        int src_index = index;
        int dest_index = index ^ (1 << cx.second);
        if(cx_dist(dest_index) > cx_dist(src_index)) {
            src_index = dest_index;
            dest_index = index;
        }
        index = dest_index;
        tracker.control_rotate_merge(src_index, dest_index);
    }
    tracker.rotate_merge(index, 0);
}

gate_sequence prepare_greedy(quantum_state state, int n_qubits) {
    StateTracker tracker = StateTracker(state);
    //phase 1
    greedy_phase1(tracker, n_qubits);
    //phase 2
    while(abs(tracker.state[0]) * abs(tracker.state[0]) < 0.95) {
        greedy_phase2_iter(tracker, n_qubits);
    }
    return tracker.gates_seq();
}
