/// @file optimizer.cxx
/// CMA-ES classes for nonlinear function optimization.
/// Copyright (c) 2017 Ralf Quast
///
/// Permission is hereby granted, free of charge, to any person obtaining a copy
/// of this software and associated documentation files (the "Software"), to deal
/// in the Software without restriction, including without limitation the rights
/// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
/// copies of the Software, and to permit persons to whom the Software is
/// furnished to do so, subject to the following conditions:
///
/// The above copyright notice and this permission notice shall be included in all
/// copies or substantial portions of the Software.
///
/// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
/// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
/// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
/// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
/// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
/// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
/// SOFTWARE.
#include <cmath>

#include "optimizer.h"

especia::Optimizer::Builder::Builder() : weights(parent_number) {
    set_strategy_parameters();
}

especia::Optimizer::Builder::~Builder() {

}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_problem_dimension(N_type n) {
    this->n = n;
    set_strategy_parameters();
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_parent_number(N_type parent_number) {
    this->parent_number = parent_number;
    set_strategy_parameters();
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_population_size(N_type population_size) {
    this->population_size = population_size;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_covariance_update_modulus(N_type update_modulus) {
    this->update_modulus = update_modulus;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_accuracy_goal(R_type accuracy_goal) {
    this->accuracy_goal = accuracy_goal;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_random_seed(W_type seed) {
    this->random_seed = seed;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_stop_generation(L_type stop_generation) {
    this->stop_generation = stop_generation;
    return *this;
}

especia::Optimizer especia::Optimizer::Builder::build() {
    return Optimizer(*this);
}

void especia::Optimizer::Builder::set_strategy_parameters() {
    using std::log;
    using std::max;
    using std::min;
    using std::sqrt;

    weights.resize(parent_number);

    for (N_type i = 0; i < parent_number; ++i) {
        weights[i] = log((parent_number + 1.0) / (i + 1));
    }

    wv = sq(weights.sum()) / weights.apply(sq).sum();
    cs = (wv + 2.0) / (wv + n + 3.0);
    cc = 4.0 / (n + 4.0);

    acov = 1.0 / wv;
    ccov = acov * (2.0 / sq(n + sqrt(2.0))) + (1.0 - acov) * min<R_type>(1.0, (2.0 * wv - 1.0) / (sq(n + 2.0) + wv));

    step_size_damping = cs + 1.0 + 2.0 * max<R_type>(0.0, sqrt((wv - 1.0) / (n + 1.0)) - 1.0);
}

especia::Optimizer::Result::Result(N_type n,
                                   const valarray<R_type> &x_in,
                                   const valarray<R_type> &d_in,
                                   R_type s_in)
        : x(x_in), d(d_in), s(s_in), z(0.0, n), B(0.0, sq(n)), C(0.0, sq(n)), pc(0.0, n), ps(0.0, n) {
    for (N_type i = 0, ii = 0; i < n; ++i, ii += n + 1) {
        B[ii] = 1.0;
        C[ii] = d[i] * d[i];
    }

    y = 0.0;

    optimized = false;
    underflow = false;

    g = 0;
}

especia::Optimizer::Result::~Result() {

}

especia::Optimizer::Optimizer(const especia::Optimizer::Builder &builder)
        : config(builder), decompose(builder.get_problem_dimension()), deviate(builder.get_random_seed()) {

}

especia::Optimizer::~Optimizer() {

}
