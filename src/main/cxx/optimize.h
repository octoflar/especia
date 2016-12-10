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
     * @tparam number The number type.
     * @tparam comparator The direct number comparator type.
     */
    template<class number, class comparator>
    class indirect_comparator {
    public:
        /**
         * Constructs a new indirect comparator to compare a set of fitness values.
         *
         * @param[in] f The fitness values.
         * @param[in] c The direct number comparator.
         */
        indirect_comparator(const std::valarray<number> &f, const comparator &c) : fitness(f), comp(c) {
        }

        /**
         * Destructor.
         */
        ~indirect_comparator() {
        }

        /**
         * The indirect comparation.
         *
         * @param[in] i An index into the set of fitness values.
         * @param[in] j An index into the set of fitness values.
         * @return the direct comparator result.
         */
        bool operator()(const size_t &i, const size_t &j) const {
            return comp(fitness[i], fitness[j]);
        }

    private:
        const std::valarray<number> &fitness;
        const comparator &comp;
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
     * @tparam function  The function type.
     * @tparam validator The validator type.
     * @tparam generator The strategy to generate random normal deviates.
     * @tparam decomposer The strategy to perform the symmetric eigenvalue decomposition.
     * @tparam comparator The strategy to compare fitness values.
     * @param[in] f The objective function.
     * @param[in] constraint The prior constraints.
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
     * @param[out] yw The value of the objective function at @c xw.
     * @param[out] optimized Set to @c true when the optimization has converged.
     * @param[out] underflow Set to @c true when the mutation variance is too small.
     * @param[in] random The random number generator.
     * @param[in] decomp The eigenvalue decomposition.
     * @param[in] comp The comparator to compare fitness values.
     */
    template<class function, class validator, class generator, class decomposer, class comparator>
    void optimize(const function &f,
                  const validator &constraint,
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
                  generator &random, decomposer &decomp, const comparator &comp) {
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
                        const double z = random();

                        for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                            u[k][i] = uw[i] + z * BD[ij];
                            v[k][i] = vw[i] + z * B[ij];
                            x[k][i] = xw[i] + u[k][i] * step_size; // Hansen and Ostermeier (2001), Eq. (13)
                        }
                    } while (constraint.reject(&x[k][0], n));
                    uw = u[k];
                    vw = v[k];
                }
            }
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t k = 0; k < population_size; ++k) {
                fitness[k] = f(&x[k][0], n);
                index[k] = k;
            }
            partial_sort(&index[0], &index[parent_number], &index[population_size],
                         indirect_comparator<double, comparator>(fitness, comp));
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
                decomp(C, B, d, n);

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
            if (optimized)
                break;
        }

        yw = f(xw, n);
    }


    /**
     * Function to scale the global step size to compute standard uncertainties and covariance.
     *
     * Compute the standard covariance along the major principal axis from the curvature of a
     * parabola through three point around the minimum.
     *
     * @tparam function The function type.
     * @param[in] f The objective function.
     * @param[in] x The parameter values.
     * @param[in] n The number of parameter values.
     * @param[in] d The local step sizes
     * @param[in] B The rotation matrix.
     * @param[in,out] s The global step size.
     */
    template<class function>
    void scale_step_size(const function &f, const double x[], size_t n, const double d[], const double B[], double &s) {
        using std::abs;
        using std::valarray;

        const double a = 100.0 * s;

        valarray<double> p(x, n);
        valarray<double> q(x, n);
        for (size_t i = 0, j = n - 1, ij = j; i < n; ++i, ij += n) {
            p[i] += a * B[ij] * d[j];
            q[i] += a * B[ij] * d[j];
        }

        const double zx = f(&x[0], n);
        const double zp = f(&p[0], n);
        const double zq = f(&q[0], n);
        // compute the covariance along the major principal axis by means of a parabola
        s = a / sqrt(abs(2.0 * (zp - zx) - (zp - zq)));
    }


    /**
     * A strict bound constraint.
     *
     * @tparam number The number type.
     */
    template<class number>
    class bound_constraint {
    public:
        /**
         * Constructs a new bound constraint.
         *
         * @param lower_bounds The lower bounds.
         * @param upper_bounds The upper bounds.
         * @param n The number of bounds.
         */
        bound_constraint(const number lower_bounds[], const number upper_bounds[], size_t n)
                : a(lower_bounds, n), b(upper_bounds, n) {
        }

        /**
         * Destructor.
         */
        ~bound_constraint() {
        }

        /**
         * Tests if a given parameter vector violates the bound constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to test.
         * @return @c true, if the parameter vector violates the bound constraint.
         */
        bool reject(const number x[], size_t n) const {
            for (size_t i = 0; i < n; ++i) {
                if (x[i] < a[i] || x[i] > b[i]) {
                    return true;
                }
            }
            return false;
        }

    private:
        const std::valarray<number> a;
        const std::valarray<number> b;
    };


    /**
     * No constraint.
     *
     * @tparam number The number type.
     */
    template<class number>
    class no_constraint {
    public:
        /**
         * Constructor.
         */
        no_constraint() {
        }

        /**
         * Destructor.
         */
        ~no_constraint() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param x The parameter vector.
         * @param n The number of parameters to test.
         * @return always @c false.
         */
        bool reject(const number x[], size_t n) const {
            return false;
        }
    };

}

#endif // ESPECIA_OPTIMIZE_H
