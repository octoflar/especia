/// @file optimize.h
/// CMA-ES function templates for nonlinear function optimization.
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
#ifndef ESPECIA_OPTIMIZE_H
#define ESPECIA_OPTIMIZE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>
#include <valarray>
#include <vector>

#include "base.h"

namespace especia {

    /**
     * An indirect comparing of indexes.
     *
     * @tparam T The base value type.
     * @tparam Compare The strategy to compare base values directly.
     */
    template<class T, class Compare>
    class Index_Compare {
    public:
        /**
         * Constructs a new index comparing.
         *
         * @param[in] v The base values.
         * @param[in] c The direct base value comparing.
         */
        Index_Compare(const std::valarray<T> &v, const Compare &c)
                : values(v), compare(c) {
        }

        /**
         * The destructor.
         */
        ~Index_Compare() {
        }

        /**
         * The index comparing operator.
         *
         * @param[in] i An index into the set of base values.
         * @param[in] j An index into the set of base values.
         * @return the result of comparing the indexed base values directly.
         */
        bool operator()(const natural &i, const natural &j) const {
            return compare(values[i], values[j]);
        }

    private:
        const std::valarray<T> &values;
        const Compare &compare;
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
     * @tparam Compare The strategy to compare fitness.
     * @tparam Tracing The tracer type.
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
     * @param[out] yw The fitness at @c xw.
     * @param[out] optimized Set to @c true when the optimization has converged.
     * @param[out] underflow Set to @c true when the mutation variance is too small.
     * @param[in] deviate The random number generator.
     * @param[in] decompose The eigenvalue decomposition.
     * @param[in] compare The comparator to compare fitness.
     * @param[in] tracer The tracer.
     */
    template<class F, class Constraint, class Deviate, class Decompose, class Compare, class Tracing>
    void optimize(const F &f,
                  const Constraint &constraint,
                  natural n,
                  natural parent_number,
                  natural population_size,
                  const real w[],
                  real step_size_damping,
                  real cs,
                  real cc,
                  real ccov,
                  real acov,
                  natural update_modulus,
                  real accuracy_goal,
                  natural stop_generation,
                  natural &g,
                  real xw[],
                  real &step_size,
                  real d[],
                  real B[],
                  real C[],
                  real ps[],
                  real pc[],
                  real &yw,
                  bool &optimized,
                  bool &underflow,
                  const Deviate &deviate, const Decompose &decompose, const Compare &compare, const Tracing &tracer) {
        using std::accumulate;
        using std::exp;
        using std::inner_product;
        using std::numeric_limits;
        using std::partial_sort;
        using std::sqrt;
        using std::thread;
        using std::valarray;
        using std::vector;

        const real expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(real(n));
        const real max_covariance_matrix_condition = 0.01 / numeric_limits<real>::epsilon();
        const real csu = sqrt(cs * (2.0 - cs));
        const real ccu = sqrt(cc * (2.0 - cc));
        const real ws = accumulate(w, w + parent_number, 0.0);
        const real cw = ws / sqrt(inner_product(w, w + parent_number, w, 0.0));

        valarray<real> uw(n);
        valarray<real> vw(n);
        valarray<valarray<real>> u(uw, population_size);
        valarray<valarray<real>> v = u;
        valarray<valarray<real>> x = u;

        valarray<real> y(population_size);
        valarray<natural> indexes(population_size);

        valarray<real> BD(B, n * n);

        while (g < stop_generation) {
            for (natural j = 0; j < n; ++j) {
                for (natural i = 0, ij = j; i < n; ++i, ij += n) {
                    BD[ij] = B[ij] * d[j];
                }
            }

            // Generate a new population of object parameter vectors,
            // sorted indirectly by fitness
            for (natural k = 0; k < population_size; ++k) {
                uw = 0.0;
                vw = 0.0;
                for (natural j = 0; j < n; ++j) {
                    do {
                        const real z = deviate();

                        for (natural i = 0, ij = j; i < n; ++i, ij += n) {
                            u[k][i] = uw[i] + z * BD[ij];
                            v[k][i] = vw[i] + z * B[ij];
                            x[k][i] = xw[i] + u[k][i] * step_size; // Hansen and Ostermeier (2001), Eq. (13)
                        }
                    } while (constraint.is_violated(&x[k][0], n));
                    uw = u[k];
                    vw = v[k];
                }
            }
#ifdef _OPENMP
#pragma omp parallel for
            for (natural k = 0; k < population_size; ++k) {
                y[k] = f(&x[k][0], n) + constraint.cost(&x[k][0], n);
                indexes[k] = k;
            }
#else // C++-11
            vector<thread> threads; threads.reserve(population_size);
            for (natural k = 0; k < population_size; ++k) {
                threads.push_back(
                        thread([k, &f, &constraint, &x, n, &y]() {
                            y[k] = f(&x[k][0], n) + constraint.cost(&x[k][0], n);
                        })
                );
            }
            for (natural k = 0; k < population_size; ++k) {
                threads[k].join();
                indexes[k] = k;
            }
#endif
            partial_sort(&indexes[0], &indexes[parent_number], &indexes[population_size],
                         Index_Compare<real, Compare>(y, compare));
            ++g;

            // Check the mutation variance
            underflow = (y[indexes[0]] == y[indexes[parent_number]]);
            if (!underflow)
                for (natural i = 0, ij = g % n; i < n; ++i, ij += n) {
                    underflow = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                    if (!underflow) {
                        break;
                    }
                }
            if (underflow) {
                break;
            }

            // Recombine the best individuals
            for (natural i = 0; i < n; ++i) {
                uw[i] = vw[i] = xw[i] = 0.0;
                for (natural j = 0; j < parent_number; ++j) {
                    uw[i] += w[j] * u[indexes[j]][i];
                    vw[i] += w[j] * v[indexes[j]][i];
                    xw[i] += w[j] * x[indexes[j]][i];
                }
                uw[i] /= ws;
                vw[i] /= ws;
                xw[i] /= ws;
            }

            real s = 0.0;
            // Adapt the covariance matrix and the step size according to Hansen and Ostermeier (2001)
            // and Hansen et al. (2003)
            for (natural i = 0, i0 = 0; i < n; ++i, i0 += n) {
                pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i]; // ibd. (2001), Eq. (14)
                if (ccov > 0.0) {
                    // BD is not used anymore and can be overwritten
                    valarray<real> &Z = BD;

                    for (natural j = 0, ij = i0; j <= i; ++j, ++ij) {
                        Z[ij] = 0.0;
                        for (natural k = 0; k < parent_number; ++k) {
                            Z[ij] += w[k] * (u[indexes[k]][i] * u[indexes[k]][j]);
                        }
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
                decompose(C, B, d);

                real t;
                // Adjust the condition of the covariance matrix and recompute the
                // local step sizes
                if ((t = d[n - 1] / max_covariance_matrix_condition - d[0]) > 0.0) {
                    for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                        C[ii] += t;
                        d[i] += t;
                    }
                }
                for (natural i = 0; i < n; ++i) {
                    d[i] = sqrt(d[i]);
                }
            }

            // Check if the optimization is completed
            for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                optimized = (sq(step_size) * C[ii] <
                        sq(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
                if (!optimized) {
                    break;
                }
            }
            if (optimized or tracer.is_enabled(g)) {
                tracer.trace(g, f(xw, n) + constraint.cost(xw, n), step_size * d[0], step_size * d[n - 1]);
            }
            if (optimized) {
                break;
            }
        }

        yw = f(xw, n) + constraint.cost(xw, n);
    }

    /**
     * Yields the paramater standard uncertainties.
     *
     * Computes the standard variance along the minor principal axis from the curvature of a
     * parabola through three points around the minimum. The global step size is rescaled to
     * approximate the standard covariance matrix by the product of the squared global step
     * size with the optimized covariance matrix.
     *
     * @tparam F The function type.
     * @tparam Constraint The constraint type.
     *
     * @param[in] f The objective function.
     * @param[in] constraint The constraint on parameter values.
     * @param[in] n The number of parameter values.
     * @param[in] x The parameter values.
     * @param[in] d The local step sizes
     * @param[in] B The rotation matrix.
     * @param[in] B The covariance matrix.
     * @param[in] s The global step size.
     * @param[out] z The parameter uncertainties.
     */
    template<class F, class Constraint>
    void postopti(const F &f, const Constraint &constraint, natural n,
                  const real x[],
                  const real d[],
                  const real B[],
                  const real C[],
                  real s,
                  real z[]) {
        using std::abs;
        using std::sqrt;
        using std::valarray;
        using std::thread;

        const real zx = f(&x[0], n) + constraint.cost(&x[0], n);

        real a = 0.0;
        real b = 0.0;
        real c = s;

        do {
            // Compute two steps along the line of least variance in opposite directions
            valarray<real> p(x, n);
            valarray<real> q(x, n);
            for (natural i = 0, j = 0, ij = j; i < n; ++i, ij += n) {
                p[i] += c * B[ij] * d[j];
                q[i] -= c * B[ij] * d[j];
            }
            real zp;
            real zq;
#ifdef _OPENMP
#pragma omp parallel
            {
#pragma omp sections
                {
#pragma omp section
                    zp = f(&p[0], n) + constraint.cost(&p[0], n);
#pragma omp section
                    zq = f(&q[0], n) + constraint.cost(&q[0], n);
                }
            }
#else // C++-11
            thread tp([&f, &constraint, &p, n, &zp]() { zp = f(&p[0], n) + constraint.cost(&p[0], n); });
            thread tq([&f, &constraint, &q, n, &zq]() { zq = f(&q[0], n) + constraint.cost(&q[0], n); });
            tp.join();
            tq.join();
#endif
            // Compute the rescaled global step size
            s = c / sqrt(abs((zp + zq) - (zx + zx)));

            // Make a smaller or larger computation step in the next iteration
            if (abs(zq - zx) < 0.1) {
                a = c;
                c = c * 1.618;
            } else {
                b = c;
                c = c * 0.618;
            }
        } while (a == 0.0 or b == 0.0); // the computation step is too small or too large

        for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
            z[i] = s * sqrt(C[ii]);
        }
    }

}

#endif // ESPECIA_OPTIMIZE_H
