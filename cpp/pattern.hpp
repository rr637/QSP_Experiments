#pragma once
#include "simulator.hpp"
#include <string>

//middle_start is inclusive, middle_end is not inclusive. So the binary string
//pattern "1101**011" would be represented as:
//    lower -> 3
//    upper -> 13
//    middle_start -> 3
//    middle_end -> 5
//
class Pattern {
    public:
        Pattern(int lower, int upper, int middle_start, int middle_end);
        int lower;
        int upper;
        int middle_start;
        int middle_end;
        int trit_at(int index) const;
        Pattern bit_flip(int index) const;
        std::string to_string(int n_qubits) const;
};

quantum_state substate(const quantum_state& state, const Pattern& pattern);
