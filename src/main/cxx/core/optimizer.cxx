//! @file optimizer.cxx
//! CMA-ES classes for nonlinear function optimization.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <cmath>

#include "optimizer.h"

especia::Optimizer::Builder::Builder() : weights(parent_number) {
    with_strategy_parameters();
}

especia::Optimizer::Builder::~Builder() = default;

especia::Optimizer::Builder &especia::Optimizer::Builder::with_defaults() {
    return with_problem_dimension().
            with_parent_number().
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
    const natural population_size = 2 * parent_number;

    if (this->parent_number != parent_number) {
        this->parent_number = parent_number;
        with_strategy_parameters();
    }
    if (this->population_size != population_size) {
        with_population_size(population_size);
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

especia::Optimizer::Builder &especia::Optimizer::Builder::with_random_seed(word64 seed) {
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

especia::Optimizer::Result::~Result() = default;

especia::Optimizer::Optimizer(const especia::Optimizer::Builder &builder)
        : config(builder), decompose(builder.get_problem_dimension()), deviate(builder.get_random_seed()) {

}

especia::Optimizer::~Optimizer() = default;
