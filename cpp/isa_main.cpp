
#include "greedy3.hpp"
#include "simulator.hpp"
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
  if(abs(norm(state) - 1) > 0.001) {
    throw std::runtime_error("Non_normalized input state!")
    ;
  }
  gate_sequence gates = greedy3_prepare(state, n_qubits, fidelity);
  std::ofstream output(argv[1]);
  for(const Gate& gate : gates) {
    output << gate.to_string() << std::endl;
  }
  return 0;
}
