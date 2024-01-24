#include "pattern.hpp"
#include "simulator.hpp"
#include <iostream>

//int main() {
//    quantum_state state = random_state(4);
//    print_state(state, 4);
//    Pattern p = Pattern(0, 1, 1, 3);
//    quantum_state subs = substate(state, p);
//    print_state(subs, 4);
//}

Pattern::Pattern(int lower, int upper, int middle_start, int middle_end) : 
    lower(lower), upper(upper), middle_start(middle_start), middle_end(middle_end) {
    if(middle_end < middle_start) throw std::runtime_error("???");
}

int Pattern::trit_at(int index) const {
    if(index < middle_start) return (this->lower >> index) % 2;
    if(index >= middle_end) return (this->upper >> (index - middle_end)) % 2;
    return -1;
}

Pattern Pattern::bit_flip(int index) const {
    if(index < middle_start) {
        int new_lower = this->lower ^ (1 << index);
        return Pattern(new_lower, this->upper, this->middle_start, 
          this->middle_end);
    }
    if(index >= middle_end) {
        int new_upper = this->upper ^ (1 << (index - middle_end));
        return Pattern(this->lower, new_upper, this->middle_start,
          this->middle_end);
    }
    throw std::runtime_error("cannot flip a wild bit");
}

std::string Pattern::to_string(int n_qubits) const {
    std::string output = "";
    for(int i = 0; i < this->middle_start; i++) {
        output = std::to_string((this->lower >> i) & 1) + output;
    }
    for(int i = this->middle_start; i < this->middle_end; i++) {
        output = "*" + output;
    }
    for(int i = this->middle_end; i < n_qubits; i++) {
        output = std::to_string((this->upper >> (i - this->middle_end)) & 1) 
          + output;
    }
    return output;
}
//Returns non-normalized output
quantum_state substate(const quantum_state& state, const Pattern& pattern) {
    quantum_state output;
    int base = (pattern.upper << pattern.middle_end) + pattern.lower;
    int term = 1 << (pattern.middle_end - pattern.middle_start);
    for(int i = 0; i < 1 << (pattern.middle_end - pattern.middle_start); i++) {
        output.push_back(state[base + (i << pattern.middle_start)]);
    }
    return output;
}

