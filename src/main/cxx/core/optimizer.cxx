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
    with_strategy_parameters();
}

especia::Optimizer::Builder::~Builder() {

}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_defaults() {
    return with_problem_dimension().
            with_covariance_update_modulus().
            with_accuracy_goal().
            with_stop_generation().
            with_random_seed();
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_problem_dimension(natural n) {
    if (n != this->n) {
        this->n = n;
        with_strategy_parameters();
    }
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_parent_number(natural parent_number) {
    if (parent_number != this->parent_number) {
        this->parent_number = parent_number;
        with_population_size(2 * parent_number);
        with_strategy_parameters();
    }
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_population_size(natural population_size) {
    this->population_size = population_size;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_covariance_update_modulus(natural update_modulus) {
    this->update_modulus = update_modulus;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_accuracy_goal(real accuracy_goal) {
    this->accuracy_goal = accuracy_goal;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_random_seed(word seed) {
    this->random_seed = seed;
    return *this;
}

especia::Optimizer::Builder &especia::Optimizer::Builder::with_stop_generation(natural stop_generation) {
    this->stop_generation = stop_generation;
    return *this;
}

especia::Optimizer especia::Optimizer::Builder::build() {
    return Optimizer(*this);
}

void especia::Optimizer::Builder::with_strategy_parameters() {
    using std::log;
    using std::max;
    using std::min;
    using std::sqrt;

    weights.resize(parent_number);

    for (natural i = 0; i < parent_number; ++i) {
        weights[i] = log((parent_number + 0.5) / (i + 1));
    }

    wv = sq(weights.sum()) / weights.apply(sq).sum();
    cs = (2.0 + wv) / (5.0 + n + wv);
    cc = (4.0 + wv / n) / (4.0 + n + 2.0 * wv / n);

    acov = 2.0 / (sq(n + 1.3) + wv);
    ccov = min<real>(1.0 - acov, 2.0 * (wv - 2.0 + 1.0 / wv) / (sq(n + 2.0) + wv));
    step_size_damping = cs + 1.0 + 2.0 * max<real>(0.0, sqrt((wv - 1.0) / (n + 1.0)) - 1.0);
}

especia::Optimizer::Result::Result(natural n,
                                   const valarray<real> &x_in,
                                   const valarray<real> &d_in,
                                   real s_in)
        : x(x_in), d(d_in), s(s_in), z(0.0, n), B(0.0, sq(n)), C(0.0, sq(n)), pc(0.0, n), ps(0.0, n) {
    for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
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
