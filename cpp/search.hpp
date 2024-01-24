
#include <vector>
#include "simulator.hpp"

std::vector<gate_sequence> prepare_search(quantum_state state, int n_qubits, int max_cx, double target_fidelity);
