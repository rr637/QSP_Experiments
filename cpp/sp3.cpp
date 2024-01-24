#include "sp3.hpp"
#include "simulator.hpp"
#include "state_tracker.hpp"
#include <vector>
#include <array>
#include <complex>
#include <iostream>

void sp3(StateTracker& tracker, int q0, int q1, int q2);

//int main() {
//    quantum_state state = random_state(5);
//    StateTracker tracker = StateTracker(state);
//    print_state(tracker.state, 5);
//    sp3(tracker, 1, 2, 3);
//    print_state(tracker.state, 5);
//}

//helper functions for handling 2x2 matrices
using mat22 = std::array<std::complex<double>, 4>;

mat22 operator* (std::complex<double> c, mat22 m) {
    mat22 output {m[0] * c, m[1] * c, m[2] * c, m[3] * c};
    return output;
}

mat22 operator* (mat22 m, std::complex<double> c) {
    return c * m;
}

mat22 operator* (mat22 m1, mat22 m2) {
    mat22 output {
        m1[0] * m2[0] + m1[1] * m2[2],
        m1[0] * m2[1] + m1[1] * m2[3],
        m1[2] * m2[0] + m1[3] * m2[2],
        m1[2] * m2[1] + m1[3] * m2[3]
    };
    return output;
}

mat22 conj_transpose(mat22 m) {
    mat22 output {conj(m[0]), conj(m[2]), conj(m[1]), conj(m[3])};
    return output;
}

mat22 phase_mat(double phase) {
    mat22 output {1, 0, 0, std::complex<double>(cos(phase), sin(phase))};
    return output;
}

std::array<std::complex<double>, 2> normalize(std::array<std::complex<double>, 2> s) {
    double norm = sqrt(abs(s[0]) * abs(s[0]) + abs(s[1]) * abs(s[1]));
    std::array<std::complex<double>, 2> output {s[0] / norm, s[1] / norm};
    return output;
}

mat22 eig_x(mat22 m) {
    std::array<std::complex<double>, 2> v1 {-m[1], m[0] - 1.0};
    std::array<std::complex<double>, 2> v2 {-m[1], m[0] + 1.0};
    v1 = normalize(v1);
    v2 = normalize(v2);
    mat22 p {v1[0], v2[0], v1[1], v2[1]};
    double a = 1 / sqrt(2);
    mat22 h {a, a, a, -a};
    return p * h;
}

class SubstateView {
    public:
        SubstateView(StateTracker& tracker, int q0, int q1, int q2) :
            tracker(tracker), q0(q0), q1(q1), q2(q2) {}
        std::complex<double> operator [](int i) {
            if(i < 0 || i > 7) {
                throw new std::runtime_error("Out of bounds");
            }
            int index = 0;
            index += (i & 1) << this->q0;
            index += (i & 2) >> 1 << this->q1;
            index += (i & 4) >> 2 << this->q2;
            return tracker.state[index];
        }
    private:
        StateTracker& tracker;
        int q0, q1, q2;
};

void zyz(StateTracker& tracker, const mat22& mat, int index) {
    double global_phase = arg(mat[0]);
    double z1 = arg(mat[2]) - global_phase;
    double y = 2 * atan(abs(mat[2]) / abs(mat[0]));
    double z2 = arg(mat[1]) - global_phase + M_PI;
    tracker.apply_rz(index, z2);
    tracker.apply_ry(index, y);
    tracker.apply_rz(index, z1);
}

std::complex<double> safe_div(std::complex<double> a, std::complex<double> b) {
    if(abs(b) > 1e-10) return a / b;
    return (a + 1e-9) / (b + 1e-9);
}

std::complex<double> safe_option_div(std::complex<double> a, 
  std::complex<double> b, std::complex<double> c, std::complex<double> d) {
    if(abs(a) + abs(b) > abs(c) + abs(d)) return safe_div(a, b);
    return safe_div(c, d);
}

//zeroes out q2
void sp3_1(StateTracker& tracker, int q0, int q1, int q2) {
    SubstateView state = SubstateView(tracker, q0, q1, q2);
    std::complex<double> r;
    std::complex<double> l;
    if(abs(state[7] * state[1] - state[5] * state[3]) < 1e-9) {
        r = -conj(safe_option_div(state[3], state[1], state[7], state[5]));
        auto c = safe_option_div(state[5], state[1], state[7], state[3]);
        l = safe_div(state[4] - c * state[0], state[6] - c * state[2]);
    } else if(abs(state[6] * state[0] - state[4] * state[2]) < 1e-9) {
        l = -conj(safe_option_div(state[2], state[0], state[6], state[4]));
        auto c = safe_option_div(state[4], state[0], state[6], state[2]);
        r = safe_div(state[5] - c * state[1], state[7] - c * state[3]);
    } else {
        std::complex<double> A = state[4] * state[1] - state[0] * state[5];
        std::complex<double> B = state[4] * state[3] - state[0] * state[7];
        std::complex<double> C = state[6] * state[1] - state[2] * state[5];
        std::complex<double> D = state[6] * state[3] - state[2] * state[7];
        std::complex<double> a = -(conj(C) * D + conj(A) * B);
        std::complex<double> b = conj(A) * A - conj(B) * B + conj(C) * C - conj(D) * D;
        std::complex<double> c = C * conj(D) + A * conj(B);
        r = safe_div(-b + sqrt(b * b - 4.0 * a * c), 2.0 * a);
        l = safe_div(-(conj(C) * r + conj(D)), conj(A) * r + conj(B));
    }
    double tr = atan(abs(r));
    std::complex<double> pr = safe_div(r, abs(r));
    double tl = atan(abs(l));
    std::complex<double> pl = safe_div(l, abs(l));

    mat22 ml {cos(tl), -sin(tl) * pl, sin(tl), cos(tl) * pl};
    mat22 mr {cos(tr), -sin(tr) * pr, sin(tr), cos(tr) * pr};
    mat22 abt = ml * conj_transpose(mr);
    double global_phase = arg(abt[0]);
    double gamma = M_PI - arg(abt[3]) + global_phase;
    
    ml = phase_mat(gamma) * ml 
      * std::complex<double>(cos(global_phase), -sin(global_phase));
    abt = ml * conj_transpose(mr);
    mat22 u = eig_x(abt);
    mat22 v = conj_transpose(u) * ml;
    zyz(tracker, v, q1);
    tracker.apply_cx(q0, q1);
    zyz(tracker, u, q1);
    if(abs(state[4]) + abs(state[0]) > abs(state[5]) + abs(state[1])) {
        tracker.rotate_merge((1 << q2), 0);
    } else {
        tracker.rotate_merge((1 << q2) + (1 << q0), (1 << q0));
    }
    if(abs(state[6]) + abs(state[2]) > abs(state[7]) + abs(state[3])) {
        tracker.control_rotate_merge((1 << q2) + (1 << q1), (1 << q1));
    } else {
        tracker.control_rotate_merge((1 << q2) + (1 << q1) + (1 << q0), \
          (1 << q1) + (1 << q0));
    }
    if(abs(state[2]) + abs(state[0]) > abs(state[3]) + abs(state[1])) {
        tracker.rotate_merge((1 << q1), 0);
    } else {
        tracker.rotate_merge((1 << q1) + (1 << q0), (1 << q0));
    }
}

void sp2(StateTracker& tracker, int q0, int q1) {
//    std::cout << "Called sp2" << std::endl;
//    std::cout << "q0" << q0 << std::endl;
//    std::cout << "q1" << q1 << std::endl;
    tracker.rotate_merge((1 << q1), 0);
    tracker.control_rotate_merge((1 << q0) + (1 << q1), (1 << q0));
    tracker.rotate_merge((1 << q0), 0);
}

void sp3(StateTracker& tracker, int q0, int q1, int q2) {
    sp3_1(tracker, q0, q1, q2);
    sp2(tracker, q0, q1);
}
