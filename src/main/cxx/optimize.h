// CMA-ES function templates for nonlinear function optimization
// Copyright (c) 2016 Ralf Quast
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef ESPECIA_OPTIMIZE_H
#define ESPECIA_OPTIMIZE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <valarray>

#include "base.h"

namespace especia {

    /**
     * An indirect comparator to compare a set of fitness values.
     *
     * @tparam T The number type.
     * @tparam Compare The strategy to compare numbers directly.
     */
    template<class T, class Compare>
    class Indirect_Compare {
    public:
        /**
         * Constructs a new indirect comparator to compare a set of fitness values.
         *
         * @param[in] f The fitness values.
         * @param[in] c The direct number comparator.
         */
        Indirect_Compare(const std::valarray<T> &f, const Compare &c)
                : fitness(f), compare(c) {
        }

        /**
         * Destructor.
         */
        ~Indirect_Compare() {
        }

        /**
         * The indirect comparing.
         *
         * @param[in] i An index into the set of fitness values.
         * @param[in] j An index into the set of fitness values.
         * @return the direct comparator result.
         */
        bool operator()(const size_t &i, const size_t &j) const {
            return compare(fitness[i], fitness[j]);
        }

    private:
        const std::valarray<T> &fitness;
        const Compare &compare;
    };


    /**
     * A bound constraint.
     *
     * @tparam T The number type.
     */
    template<class T>
    class Bound_Constraint {
    public:
        /**
         * Constructs a new strict-bound prior constraint.
         *
         * @param lower_bounds The lower bounds.
         * @param upper_bounds The upper bounds.
         * @param n The number of bounds.
         */
        Bound_Constraint(const T lower_bounds[], const T upper_bounds[], size_t n)
                : a(lower_bounds, n), b(upper_bounds, n) {
        }

        /**
         * Destructor.
         */
        ~Bound_Constraint() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to test.
         * @return @c true, if the parameter vector violates the constraint.
         */
        bool test(const T x[], size_t n) const {
            for (size_t i = 0; i < n; ++i) {
                if (x[i] < a[i] || x[i] > b[i]) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Computes the cost associated with the constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to take account of.
         * @return always zero.
         */
        T cost(const T x[], size_t n) const {
            return T(0);
        }

    private:
        const std::valarray<T> a;
        const std::valarray<T> b;
    };


    /**
     * No constraint.
     *
     * @tparam T The number type.
     */
    template<class T>
    class No_Constraint {
    public:
        /**
         * Constructor.
         */
        No_Constraint() {
        }

        /**
         * Destructor.
         */
        ~No_Constraint() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to test.
         * @return always @c false.
         */
        bool test(const T x[], size_t n) const {
            return false;
        }

        /**
         * Computes the cost associated with the constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to take account of.
         * @return always zero.
         */
        T cost(const T x[], size_t n) const {
            return T(0);
        }
    };


    /**
     * No tracer.
     *
     * @tparam T The number type.
     */
    template<class T>
    class No_Tracer {
    public:

        /**
         * Constructor.
         */
        No_Tracer() {
        }

        /**
         * Destructor.
         */
        ~No_Tracer() {
        }

        /**
         * Tests if tracing is enabled.
         *
         * @param g The generation number.
         * @return always @c false.
         */
        bool is_enabled(const unsigned long &g) const {
            return false;
        }

        /**
         * Traces some state information.
         *
         * @param g The generation number.
         * @param y The value of the objective function.
         * @param min_step The minimum step size.
         * @param max_step The maximum step size.
         */
        void trace(const unsigned long &g, const T &y, const T &min_step, const T &max_step) {
        }
    };

    template<class Deviate, class Decompose, class Compare, class Tracer>
    class Optimizer_Builder;

    template<class Deviate, class Decompose, class Compare, class Tracer>
    class Optimizer {
    public:
        ~Optimizer() {

        }

    private:
        Optimizer(Optimizer_Builder<Deviate, Decompose, Compare, Tracer> *builder) {
        }

        friend class Optimizer_Builder<Deviate, Decompose, Compare, Tracer>;
    };

    template<class Deviate, class Decompose, class Compare, class Tracer>
    class Optimizer_Builder {
    public:
        Optimizer_Builder(const Deviate &dev, const Tracer &tr, size_t dim = 1)
                : deviate(dev),
                  tracer(tr),
                  decompose(Decompose(dim)),
                  compare(Compare()),
                  n(dim) {
            set_parent_number();
            set_population_size();
            set_update_modulus();
            set_accuracy_goal();
            set_stop_generation();
        }

        ~Optimizer_Builder() {

        }

        Optimizer_Builder& set_parent_number(unsigned parent_number = 4) {
            this->parent_number = parent_number;
            return *this;
        }

        Optimizer_Builder& set_population_size(unsigned population_size = 8) {
            this->population_size = population_size;
            return *this;
        }

        Optimizer_Builder& set_update_modulus(unsigned update_modulus = 1) {
            this->update_modulus = update_modulus;
            return *this;
        }

        Optimizer_Builder& set_accuracy_goal(double accuracy_goal = 1.0E-06) {
            this->accuracy_goal = accuracy_goal;
            return *this;
        }

        Optimizer_Builder& set_stop_generation(unsigned long stop_generation = 1000) {
            this->stop_generation = stop_generation;
            return *this;
        }

        Optimizer<Deviate, Decompose, Compare, Tracer> build() {
            return Optimizer<Deviate, Decompose, Compare, Tracer>(this);
        };

        const Deviate &get_deviate() const {
            return deviate;
        }

        const Tracer &get_tracer() const {
            return tracer;
        }

        const Decompose &get_decompose() const {
            return decompose;
        }

        const Compare &get_compare() const {
            return compare;
        }

        const size_t get_n() const {
            return n;
        }

        unsigned int get_parent_number() const {
            return parent_number;
        }

        unsigned int get_population_size() const {
            return population_size;
        }

        unsigned int get_update_modulus() const {
            return update_modulus;
        }

        double get_accuracy_goal() const {
            return accuracy_goal;
        }

        unsigned long get_stop_generation() const {
            return stop_generation;
        }

    private:
        const Deviate &deviate;
        const Tracer &tracer;
        const Decompose decompose;
        const Compare compare;
        const size_t n;

        unsigned int parent_number;
        unsigned int population_size;
        unsigned int update_modulus;
        double accuracy_goal;
        unsigned long stop_generation;
    };

    /**
     * Evolution strategy with covariance matrix adaption (CMA-ES) for nonlinear
     * function optimization. Based on Hansen and Ostermeier (2001).
     *
     * Further reading:
     *
     * N. Hansen, S. D. Müller, P. Koumoutsakos (2003).
     *   *Reducing the Increasing the Time Complexity of the Derandomized Evolution
     *      Strategy with Covariance Matrix Adaption (CMA-ES).*
     *   Evolutionary Computation, 11, 1, ISSN 1063-6560.
     *
     *  N. Hansen, A. Ostermeier (2001).
     *    *Completely Derandomized Self-Adaption in Evolution Strategies.*
     *    Evolutionary Computation, 9, 159, ISSN 1063-6560.
     *
     * @tparam F The function type.
     * @tparam Constraint The constraint type.
     * @tparam Deviate The strategy to generate random normal deviates.
     * @tparam Decompose The strategy to perform the symmetric eigenvalue decomposition.
     * @tparam Compare The strategy to compare fitness values.
     * @tparam Tracer The tracer type.
     *
     * @param[in] f The model function.
     * @param[in] constraint The prior constraint on the parameter values.
     * @param[in] n The number of parameters.
     * @param[in] parent_number The number of parents per generation.
     * @param[in] population_size The number of individuals per generation. Twice the parent number, at least
     * @param[in] w The recombination weights.
     * @param[in] step_size_damping The step size damping.
     * @param[in] cs The step size cumulation rate.
     * @param[in] cc The distribution cumulation rate.
     * @param[in] ccov The covariance matrix adaption rate.
     * @param[in] acov The covariance matrix adaption mixing.
     * @param[in] update_modulus The covariance matrix update modulus.
     * @param[in] accuracy_goal The accuracy goal.
     * @param[in] stop_generation The stop generation.
     * @param[in,out] g The generation number.
     * @param[in,out] xw The parameter values.
     * @param[in,out] step_size The global step size.
     * @param[in,out] d The local step sizes.
     * @param[in,out] B The rotation matrix.
     * @param[in,out] C The covariance matrix.
     * @param[in,out] ps The step size cumulation path.
     * @param[in,out] pc The distribution cumulation path.
     * @param[out] yw The value of the objective function (plus the constraint cost) at @c xw.
     * @param[out] optimized Set to @c true when the optimization has converged.
     * @param[out] underflow Set to @c true when the mutation variance is too small.
     * @param[in] deviate The random number generator.
     * @param[in] decompose The eigenvalue decomposition.
     * @param[in] compare The comparator to compare fitness values.
     * @param[in] tracer The tracer.
     */
    template<class F, class Constraint, class Deviate, class Decompose, class Compare, class Tracer>
    void optimize(const F &f,
                  const Constraint &constraint,
                  size_t n,
                  unsigned parent_number,
                  unsigned population_size,
                  const double w[],
                  double step_size_damping,
                  double cs,
                  double cc,
                  double ccov,
                  double acov,
                  unsigned update_modulus,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  unsigned long &g,
                  double xw[],
                  double &step_size,
                  double d[],
                  double B[],
                  double C[],
                  double ps[],
                  double pc[],
                  double &yw,
                  bool &optimized,
                  bool &underflow,
                  Deviate &deviate, Decompose &decompose, const Compare &compare, Tracer &tracer) {
        using std::accumulate;
        using std::exp;
        using std::inner_product;
        using std::numeric_limits;
        using std::partial_sort;
        using std::sqrt;
        using std::valarray;

        const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
        const double max_covariance_matrix_condition = 0.01 / numeric_limits<double>::epsilon();
        const double csu = sqrt(cs * (2.0 - cs));
        const double ccu = sqrt(cc * (2.0 - cc));
        const double ws = accumulate(w, w + parent_number, 0.0);
        const double cw = ws / sqrt(inner_product(w, w + parent_number, w, 0.0));

        while (g < stop_generation) {
            valarray<double> uw(n);
            valarray<double> vw(n);
            valarray<valarray<double> > u(uw, population_size);
            valarray<valarray<double> > v = u;
            valarray<valarray<double> > x = u;

            valarray<double> fitness(population_size);
            valarray<size_t> index(population_size);

            valarray<double> BD(B, n * n);
            for (size_t j = 0; j < n; ++j)
                for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                    BD[ij] *= d[j];

            // Generate a new population of object parameter vectors,
            // sorted indirectly by fitness
            for (size_t k = 0; k < population_size; ++k) {
                uw = 0.0;
                vw = 0.0;
                for (size_t j = 0; j < n; ++j) {
                    do {
                        const double z = deviate();

                        for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                            u[k][i] = uw[i] + z * BD[ij];
                            v[k][i] = vw[i] + z * B[ij];
                            x[k][i] = xw[i] + u[k][i] * step_size; // Hansen and Ostermeier (2001), Eq. (13)
                        }
                    } while (constraint.test(&x[k][0], n));
                    uw = u[k];
                    vw = v[k];
                }
            }
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t k = 0; k < population_size; ++k) {
                fitness[k] = f(&x[k][0], n) + constraint.cost(&x[k][0], n);
                index[k] = k;
            }
            partial_sort(&index[0], &index[parent_number], &index[population_size],
                         Indirect_Compare<double, Compare>(fitness, compare));
            ++g;

            // Check the mutation variance
            underflow = (fitness[index[0]] == fitness[index[parent_number]]);
            if (!underflow)
                for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                    underflow = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                    if (!underflow)
                        break;
                }
            if (underflow)
                break;

            // Recombine the best individuals
            for (size_t i = 0; i < n; ++i) {
                uw[i] = vw[i] = xw[i] = 0.0;
                for (size_t j = 0; j < parent_number; ++j) {
                    uw[i] += w[j] * u[index[j]][i];
                    vw[i] += w[j] * v[index[j]][i];
                    xw[i] += w[j] * x[index[j]][i];
                }
                uw[i] /= ws;
                vw[i] /= ws;
                xw[i] /= ws;
            }

            double s = 0.0;
            double t = 0.0;

            // Adapt the covariance matrix and the step size according to Hansen and Ostermeier (2001)
            // and Hansen et al. (2003)
            for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
                pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i]; // ibd. (2001), Eq. (14)
                if (ccov > 0.0) {
                    // BD is not used anymore and can be overwritten
                    valarray<double> &Z = BD;

                    for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                        Z[ij] = 0.0;
                        for (size_t k = 0; k < parent_number; ++k)
                            Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                        // ibd. (2003), Eq. (11)
                        C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) + (1.0 - acov) * Z[ij] / ws);
                    }
                }
                ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i]; // ibd. (2001), Eq. (16)
                s += ps[i] * ps[i];
            }
            step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0)); // ibd. (2001), Eq. (17)

            if (ccov > 0.0 and g % update_modulus == 0) {
                // Decompose the covariance matrix and sort its eigenvalues in ascending
                // order, along with eigenvectors
                decompose(C, B, d, n);

                // Adjust the condition of the covariance matrix and recompute the
                // local step sizes
                if ((t = d[n - 1] / max_covariance_matrix_condition - d[0]) > 0.0)
                    for (size_t i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                        C[ii] += t;
                        d[i] += t;
                    }
                for (size_t i = 0; i < n; ++i)
                    d[i] = sqrt(d[i]);
            }

            // Check if the optimization is completed
            for (size_t i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                optimized = (sqr(step_size) * C[ii] <
                             sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
                if (!optimized)
                    break;
            }
            if (optimized or tracer.is_enabled(g))
                tracer.trace(g, f(xw, n) + constraint.cost(xw, n), step_size * d[0], step_size * d[n - 1]);
            if (optimized)
                break;
        }

        yw = f(xw, n) + constraint.cost(xw, n);
    }


    /**
     * Evolution strategy with covariance matrix adaption (CMA-ES) for nonlinear
     * function optimization. Based on Hansen and Ostermeier (2001).
     *
     * Further reading:
     *
     * N. Hansen, S. D. Müller, P. Koumoutsakos (2003).
     *   *Reducing the Increasing the Time Complexity of the Derandomized Evolution
     *      Strategy with Covariance Matrix Adaption (CMA-ES).*
     *   Evolutionary Computation, 11, 1, ISSN 1063-6560.
     *
     *  N. Hansen, A. Ostermeier (2001).
     *    *Completely Derandomized Self-Adaption in Evolution Strategies.*
     *    Evolutionary Computation, 9, 159, ISSN 1063-6560.
     *
     * @tparam F The function type.
     * @tparam Constraint The constraint type.
     * @tparam Deviate The strategy to generate random normal deviates.
     * @tparam Decompose The strategy to perform the symmetric eigenvalue decomposition.
     * @tparam Tracer The tracer type.
     *
     * @param[in] f The model function.
     * @param[in] constraint The prior constraint on the parameter values.
     * @param[in] n The number of parameters.
     * @param[in] parent_number The number of parents per generation.
     * @param[in] population_size The number of individuals per generation. Twice the parent number, at least
     * @param[in] update_modulus The covariance matrix update modulus.
     * @param[in] accuracy_goal The accuracy goal.
     * @param[in] stop_generation The stop generation.
     * @param[in,out] g The generation number.
     * @param[in,out] x The parameter values.
     * @param[in,out] step_size The global step size.
     * @param[in,out] d The local step sizes.
     * @param[in,out] B The rotation matrix.
     * @param[in,out] C The covariance matrix.
     * @param[out] y The value of the objective function (plus the constraint cost) at @c x.
     * @param[out] optimized Set to @c true when the optimization has converged.
     * @param[out] underflow Set to @c true when the mutation variance is too small.
     * @param[in] deviate The random number generator.
     * @param[in] decompose The eigenvalue decomposition.
     * @param[in] tracer The tracer.
     */
    template<class F, class Constraint, class Deviate, class Decompose, class Tracer>
    void minimize(const F &f,
                  const Constraint &constraint,
                  size_t n,
                  unsigned parent_number,
                  unsigned population_size,
                  unsigned update_modulus,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  unsigned long &g,
                  double x[],
                  double &step_size,
                  double d[],
                  double B[],
                  double C[],
                  double &y,
                  bool &optimized,
                  bool &underflow,
                  Deviate &deviate, Decompose &decompose, Tracer &tracer) {
        using std::less;
        using std::log;
        using std::max;
        using std::min;
        using std::sqrt;
        using std::valarray;

        valarray<double> w(1.0, parent_number);
        for (size_t i = 0; i < parent_number; ++i)
            w[i] = log((parent_number + 1.0) / (i + 1));

        const double wv = sqr(w.sum()) / w.apply(sqr).sum();
        const double cs = (wv + 2.0) / (wv + n + 3.0);
        const double cc = 4.0 / (n + 4.0);
        const double acov = 1.0 / wv;
        const double ccov = acov * (2.0 / sqr(n + sqrt(2.0))) +
                            (1.0 - acov) * min(1.0, (2.0 * wv - 1.0) / (sqr(n + 2.0) + wv));
        const double step_size_damping = cs + 1.0 + 2.0 * max(0.0, sqrt((wv - 1.0) / (n + 1.0)) - 1.0);

        valarray<double> pc(0.0, n);
        valarray<double> ps(0.0, n);

        optimize(f, constraint,
                 n,
                 parent_number,
                 population_size,
                 &w[0],
                 step_size_damping,
                 cs,
                 cc,
                 ccov,
                 acov,
                 update_modulus,
                 accuracy_goal,
                 stop_generation,
                 g,
                 x,
                 step_size,
                 d,
                 B,
                 C,
                 &ps[0],
                 &pc[0],
                 y,
                 optimized,
                 underflow,
                 deviate, decompose, less<double>(), tracer);
    }


    /**
     * Function to scale the global step size to compute standard uncertainties and covariance.
     *
     * Compute the standard covariance along the major principal axis from the curvature of a
     * parabola through three point around the minimum.
     *
     * @tparam F The function type.
     * @tparam Constraint The constraint type.
     *
     * @param[in] f The objective function.
     * @param[in] constraint The prior constraint on the parameter values.
     * @param[in] x The parameter values.
     * @param[in] n The number of parameter values.
     * @param[in] d The local step sizes
     * @param[in] B The rotation matrix.
     * @param[in,out] s The global step size.
     */
    template<class F, class Constraint>
    void scale_step_size(const F &f, const Constraint &constraint, const double x[], size_t n, const double d[],
                         const double B[], double &s) {
        using std::abs;
        using std::sqrt;
        using std::valarray;

        const double a = 100.0 * s;

        valarray<double> p(x, n);
        valarray<double> q(x, n);
        for (size_t i = 0, j = n - 1, ij = j; i < n; ++i, ij += n) {
            p[i] += a * B[ij] * d[j];
            q[i] += a * B[ij] * d[j];
        }

        const double zx = f(&x[0], n) + constraint.cost(&x[0], n);
        const double zp = f(&p[0], n) + constraint.cost(&p[0], n);
        const double zq = f(&q[0], n) + constraint.cost(&q[0], n);

        s = a / sqrt(abs(2.0 * (zp - zx) - (zp - zq)));
    }

}

#endif // ESPECIA_OPTIMIZE_H
