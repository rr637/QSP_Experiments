
#include "simulator.hpp"
#include "state_tracker.hpp"
#include <vector>

gate_sequence prepare_greedy(quantum_state state, int n_qubits);
void greedy_phase1(StateTracker& tracker, int n_qubits);
void greedy_phase2_iter(StateTracker& tracker, int n_qubits);
int max_index(std::vector<double> list);
