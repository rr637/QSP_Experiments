#pragma once

#include "simulator.hpp"
#include <vector>

//Returns fidelity for [gate_sequence] trying to prepare [quantum_state]
double fitness(const quantum_state&, const gate_sequence&);
gate_sequence gradient_descent(const gate_sequence&, const quantum_state&, double step_size);
gate_sequence gradient_descent(const gate_sequence&, const quantum_state&, const quantum_state&, double step_size);
std::vector<double> compute_gradient(const gate_sequence& gates, const quantum_state& target);
std::vector<double> compute_gradient(const gate_sequence& gates, const quantum_state& start, const quantum_state& target);
