#pragma once
#include "simulator.hpp"
#include "pattern.hpp"
#include <vector>

class StateTracker {
    public:
        StateTracker(const quantum_state& s);
        void apply(const Gate& g);
        void apply_ry(int target, double angle);
        void apply_rz(int target, double angle);
        void apply_cx(int control, int target);
        void new_block();
        void undo_gate();
        void undo_block();
        void rotate_merge(int src, int dest);
        void control_rotate_merge(int src, int dest);
        void unify_phase(int src, int dest);
        void rotate_merge_chunk(Pattern src_pattern, Pattern dest_pattern);
        void control_rotate_merge_chunk(Pattern src_pattern, Pattern dest_pattern);
        gate_sequence gates_seq();
        int cx_count();
        quantum_state state;
    private:
        std::vector<gate_sequence> gates;
};

int get_control_pattern(Pattern src, int target);
int get_target_pattern(Pattern src, Pattern dest);
int get_target(int src, int dest);
int get_control(int src, int target);
