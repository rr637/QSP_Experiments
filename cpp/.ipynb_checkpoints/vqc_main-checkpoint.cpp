
#include "simulator.hpp"
#include "vqc.hpp"
#include <complex>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>

std::vector<std::string> tokenize(const std::string& s) {
    std::vector<std::string> output;
    std::vector<char> buffer;
    for(int i = 0; i < s.length(); i++) {
        if(isspace(s[i])) {
            if(buffer.size() > 0) {
                output.push_back(std::string(buffer.begin(), buffer.end()));
                buffer.clear();
            }
        } else {
            buffer.push_back(s[i]);
        }   
    }
    if(buffer.size() > 0) {
        output.push_back(std::string(buffer.begin(), buffer.end()));
    }
    return output;
}

int main(int argc, char* argv[]) {
    if(argc < 4) {
        std::cout << "Usage: ./isa_cpp [file_output] [file_input] [n_qubits] [fidelity]" << std::endl;
        return 0;
    }
    int n_qubits = std::stoi(argv[3]);
    double fidelity = 0.95;
    if(argc > 4) {
        fidelity = std::stod(argv[4]);
    }
    std::ifstream file(argv[2]);
    std::vector<std::complex<double> > state;
    for(int i = 0; i < (1 << n_qubits); i++) {
        std::string line;
        getline(file, line);
        std::vector<std::string> tokens = tokenize(line);
        double amp = std::stod(tokens[1]);
        double phase = std::stod(tokens[2]);
        state.push_back(std::complex<double>(amp * cos(phase), amp * sin(phase)));
    }
    if(abs(norm(state) - 1) > 0.0001) {
        throw std::runtime_error("Non-normalized input state!");
    }
    gate_sequence gates;
    int num_layers = 2; //TODO: UPDATE THIS
    for(int i = 0; i < n_qubits; i++) {
        gates.push_back(Gate::RY(i, 0));
        gates.push_back(Gate::RZ(i, 0));
    }
    for(int l = 0; l < num_layers; l++) {
        for(int i = l & 1; i < n_qubits - 1; i++) {
            gates.push_back(Gate::CX(i, i + 1));
            gates.push_back(Gate::RY(i, 0));
            gates.push_back(Gate::RZ(i, 0));
            gates.push_back(Gate::RY(i + 1, 0));
            gates.push_back(Gate::RZ(i + 1, 0));
        }
    }
    int iters = 0;
    while(iters < 10000 && fitness(state, gates) < fidelity) {
        iters++;
        gates = gradient_descent(gates, state, 0.1);
    }
    std::ofstream output(argv[1]);
    for(const Gate& gate : gates) {
        output << gate.to_string() << std::endl;
    }
    return 0;
}
