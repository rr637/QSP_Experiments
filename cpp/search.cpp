
#include "simulator.hpp"
#include "vqc.hpp"
#include "search.hpp"

#include <vector>
#include <algorithm>
#include <iostream>

//int main() {
//    int n_qubits = 6;
//    int n_cx = 40;
//    double total = 0.0;
//    for(int i = 0; i < 100; i++) {
//        std::cout << i << std::endl;
//        quantum_state normal = random_state(n_qubits);
//        std::vector<gate_sequence> seqs = prepare_search(normal, n_qubits, n_cx, 0.95);
//        total += cx_count(seqs.back());
//    }
//    std::cout << total / 100 << std::endl;
//    return 0;
//}

std::vector<gate_sequence> prepare_search(quantum_state state, int n_qubits, 
    int max_cx, double target_fidelity) {
    std::vector<gate_sequence> output;
    std::vector<gate_sequence> candy;
    gate_sequence seq;
    for(int i = 0; i < n_qubits; i++) {
        seq.push_back(Gate::RY(i, 0));
        seq.push_back(Gate::RZ(i, 0));
    }
    for(int i = 0; i < 1000; i++) {
        seq = gradient_descent(seq, state, 0.1);
    }
    output.push_back(seq);
    candy.push_back(seq);
    for(int cx = 0; cx < max_cx; cx++) {
        int parents = candy.size();
        for(int i = 0; i < parents; i++) {
            for(int j = 0; j < n_qubits - 1; j++) {
                gate_sequence child;
                child.push_back(Gate::RZ(j, 0));
                child.push_back(Gate::RY(j, 0));
                child.push_back(Gate::RX(j + 1, 0));
                child.push_back(Gate::RY(j + 1, 0));
                child.push_back(Gate::CX(j, j + 1));
                child.insert(child.end(), candy[i].begin(), candy[i].end());
                for(int i = 0; i < 30; i++) {
                    child = gradient_descent(child, state, 0.1);
                }
                candy.push_back(child);
            }
        }
        std::vector<double> fidelity;
        for(int i = 0; i < candy.size(); i++) {
            double f = fitness(state, candy[i]);
            fidelity.push_back(fitness(state, candy[i]));
        }
        std::vector<double> indices;
        for(int i = 0; i < candy.size(); i++) {
            indices.push_back(i);
        }
        std::sort(indices.begin(), indices.end(), [&fidelity](int i, int j) {
            return fidelity[j] < fidelity[i];
        });
        std::vector<gate_sequence> temp;
        for(int i = 0; i < 3 && i < candy.size(); i++) {
            gate_sequence survivor = candy[indices[i]];
            for(int i = 0; i < 200; i++) {
                survivor = gradient_descent(survivor, state, 0.1);
            }
            temp.push_back(survivor);
        }
        candy = temp;
        output.push_back(candy[0]);
        if(fidelity[indices[0]] > target_fidelity) break;
    }
    return output;
}
