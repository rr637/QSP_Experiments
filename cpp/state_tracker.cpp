#include "state_tracker.hpp"
#include "pattern.hpp"
#include <vector>
#include <iostream>

//int main() {
//    StateTracker tracker(random_state(4));
//    tracker.rotate_merge(13, 5);
//    tracker.rotate_merge(5, 4);
//    tracker.rotate_merge(4, 0);
//    print_state(tracker.state, 4);
//    Pattern src = Pattern(0, 2, 0, 2);
//    Pattern dest = Pattern(0, 3, 0, 2);
//    tracker.control_rotate_merge_chunk(src, dest);
//    print_state(tracker.state, 4);
//}

//Requires: src ^ dest is a power of 2
int get_target(int src, int dest) {
    int x = src ^ dest;
    int output = 0;
    while(x > 1) {
        output += 1;
        x = x >> 1;
    }
    return output;
}

int get_control(int src, int target) {
    if((src >> (target + 1)) % 2) return target + 1;
    return target - 1;
}

StateTracker::StateTracker(const quantum_state& s) {
    this->gates = std::vector<gate_sequence>();
    this->gates.emplace_back();
    this->state = s;
}

void StateTracker::apply(const Gate& g) {
    this->gates.back().push_back(g);
    this->state = apply_gate(g, this->state);
}

void StateTracker::apply_ry(int target, double angle) {
    this->apply(Gate::RY(target, angle));
}

void StateTracker::apply_rz(int target, double angle) {
    this->apply(Gate::RZ(target, angle));
}

void StateTracker::apply_cx(int control, int target) {
    this->apply(Gate::CX(control, target));
}

void StateTracker::new_block() {
    this->gates.emplace_back();
}

void StateTracker::undo_gate() {
    gate_sequence& block = this->gates.back();
    if(block.size() == 0) {
        throw std::runtime_error("Cannot remove gate from empty block");
    }
    Gate g = block.back();
    this->state = apply_gate(g.inverse(), this->state);
    block.pop_back();
}

void StateTracker::undo_block() {
    gate_sequence& block = this->gates.back();
    quantum_state temp = this->state;
    for(int i = block.size() - 1; i >= 0; i--) {
        temp = apply_gate(block[i].inverse(), temp);
    }
    this->state = temp;
    if(this->gates.size() == 1) this->gates[0].clear();
    else this->gates.pop_back();
}

void StateTracker::rotate_merge(int src, int dest) {
    this->unify_phase(src, dest);
    int target = get_target(src, dest);
    double angle = -2 * atan(abs(this->state[src] / this->state[dest]));
    if(dest > src) angle *= -1;
    this->apply_ry(target, angle);
}

void StateTracker::control_rotate_merge(int src, int dest) {
    this->unify_phase(src, dest);
    int target = get_target(src, dest);
    int control = get_control(src, target);
    double angle = atan(abs(this->state[dest] / this->state[src]));
    if(dest > src) angle *= -1;
    this->apply_ry(target, angle);
    this->apply_cx(control, target);
    this->apply_ry(target, -angle);
}

void StateTracker::unify_phase(int src, int dest) {
    int target = get_target(src, dest);
    double angle = arg(this->state[dest]) - arg(this->state[src]);
    if(dest > src) angle *= -1;
    this->apply_rz(target, angle);
}

int get_target_pattern(Pattern src, Pattern dest) {
    if(src.middle_start != dest.middle_start 
      || src.middle_end != dest.middle_end) {
        throw std::runtime_error("mismatched patterns");
    }
    if(src.lower == dest.lower) {
        return get_target(src.upper, dest.upper) + src.middle_end;
    }
    return get_target(src.lower, dest.lower);
}

int get_control_pattern(Pattern src, int target) {
    if(src.trit_at(target + 1) == 1) return target + 1;
    return target - 1;
}

//returns phi, theta
std::pair<double, double> chunk_angle(const quantum_state& state, 
  Pattern src, Pattern dest) {
    int target = get_target_pattern(src, dest);
    bool down_rotate = src.trit_at(target) == 1;
    quantum_state v0;
    quantum_state v1;
    if(down_rotate) {
        v0 = substate(state, dest);
        v1 = substate(state, src);
    } else {
        v0 = substate(state, src);
        v1 = substate(state, dest);
    }
    double v0_norm = norm(v0);
    double v1_norm = norm(v1);
    double v0_mag = v0_norm * v0_norm;
    double v1_mag = v1_norm * v1_norm;
    std::complex<double> cross = 0;
    for(int i = 0; i < v0.size(); i++) cross += conj(v0[i]) * v1[i];
    double cross_mag = abs(cross);
    double phi = -arg(cross);
    double A = (v0_mag - v1_mag) / 2;
    double theta = -acos(A / sqrt(A * A + cross_mag * cross_mag));
    if(!down_rotate) theta += M_PI;
    return std::pair<double, double>(phi, theta);
}

void StateTracker::rotate_merge_chunk(Pattern src, Pattern dest) {
    std::pair<double, double> angles = chunk_angle(this->state, src, dest);
    int target = get_target_pattern(src, dest);
    this->apply_rz(target, angles.first);
    this->apply_ry(target, angles.second);
}

void StateTracker::control_rotate_merge_chunk(Pattern src, Pattern dest) {
    std::pair<double, double> angles = chunk_angle(this->state, src, dest);
    double theta = (angles.second - M_PI) / 2;
    int target = get_target_pattern(src, dest);
    int control = get_control_pattern(src, target);
    this->apply_rz(target, angles.first);
    this->apply_ry(target, theta);
    this->apply_cx(control, target);
    this->apply_ry(target, -theta);
}

gate_sequence StateTracker::gates_seq() {
    gate_sequence output;
    for(const gate_sequence& seq : this->gates) {
        for(const Gate& g : seq) output.push_back(g);
    }
    return output;
}

int StateTracker::cx_count() {
    int output = 0;
    for(const gate_sequence& seq : this->gates) {
        for(const Gate& g : seq) {
            if(g.is_cx()) output++;
        }
    }
    return output;
}
