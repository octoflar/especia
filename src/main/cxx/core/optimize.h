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
         * Destructor.
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
        bool operator()(const N_type &i, const N_type &j) const {
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
     * @param[out] yw The fitness at @c xw.
     * @param[out] optimized Set to @c true when the optimization has converged.
     * @param[out] underflow Set to @c true when the mutation variance is too small.
     * @param[in] deviate The random number generator.
     * @param[in] decompose The eigenvalue decomposition.
     * @param[in] compare The comparator to compare fitness.
     * @param[in] tracer The tracer.
     */
    template<class F, class Constraint, class Deviate, class Decompose, class Compare, class Tracer>
    void optimize(const F &f,
                  const Constraint &constraint,
                  N_type n,
                  N_type parent_number,
                  N_type population_size,
                  const R_type w[],
                  R_type step_size_damping,
                  R_type cs,
                  R_type cc,
                  R_type ccov,
                  R_type acov,
                  N_type update_modulus,
                  R_type accuracy_goal,
                  L_type stop_generation,
                  L_type &g,
                  R_type xw[],
                  R_type &step_size,
                  R_type d[],
                  R_type B[],
                  R_type C[],
                  R_type ps[],
                  R_type pc[],
                  R_type &yw,
                  bool &optimized,
                  bool &underflow,
                  Deviate &deviate, Decompose &decompose, const Compare &compare, Tracer &tracer) {
        using std::accumulate;
        using std::exp;
        using std::inner_product;
        using std::numeric_limits;
        using std::partial_sort;
        using std::sqrt;
        using std::thread;
        using std::valarray;
        using std::vector;

        const R_type expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(R_type(n));
        const R_type max_covariance_matrix_condition = 0.01 / numeric_limits<R_type>::epsilon();
        const R_type csu = sqrt(cs * (2.0 - cs));
        const R_type ccu = sqrt(cc * (2.0 - cc));
        const R_type ws = accumulate(w, w + parent_number, 0.0);
        const R_type cw = ws / sqrt(inner_product(w, w + parent_number, w, 0.0));

        while (g < stop_generation) {
            valarray<R_type> uw(n);
            valarray<R_type> vw(n);
            valarray<valarray<R_type> > u(uw, population_size);
            valarray<valarray<R_type> > v = u;
            valarray<valarray<R_type> > x = u;

            valarray<R_type> y(population_size);
            valarray<N_type> indexes(population_size);

            valarray<R_type> BD(B, n * n);
            for (N_type j = 0; j < n; ++j) {
                for (N_type i = 0, ij = j; i < n; ++i, ij += n) {
                    BD[ij] *= d[j];
                }
            }

            // Generate a new population of object parameter vectors,
            // sorted indirectly by fitness
            for (N_type k = 0; k < population_size; ++k) {
                uw = 0.0;
                vw = 0.0;
                for (N_type j = 0; j < n; ++j) {
                    do {
                        const R_type z = deviate();

                        for (N_type i = 0, ij = j; i < n; ++i, ij += n) {
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
            for (N_type k = 0; k < population_size; ++k) {
                y[k] = f(&x[k][0], n) + constraint.cost(&x[k][0], n);
                indexes[k] = k;
            }
#else // C++-11
            vector<thread> threads; threads.reserve(population_size);
            for (N_type k = 0; k < population_size; ++k) {
                threads.push_back(
                        thread([k, &f, &constraint, &x, n, &y]() {
                            y[k] = f(&x[k][0], n) + constraint.cost(&x[k][0], n);
                        })
                );
            }
            for (N_type k = 0; k < population_size; ++k) {
                threads[k].join();
                indexes[k] = k;
            }
#endif
            partial_sort(&indexes[0], &indexes[parent_number], &indexes[population_size],
                         Index_Compare<R_type, Compare>(y, compare));
            ++g;

            // Check the mutation variance
            underflow = (y[indexes[0]] == y[indexes[parent_number]]);
            if (!underflow)
                for (N_type i = 0, ij = g % n; i < n; ++i, ij += n) {
                    underflow = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                    if (!underflow) {
                        break;
                    }
                }
            if (underflow) {
                break;
            }

            // Recombine the best individuals
            for (N_type i = 0; i < n; ++i) {
                uw[i] = vw[i] = xw[i] = 0.0;
                for (N_type j = 0; j < parent_number; ++j) {
                    uw[i] += w[j] * u[indexes[j]][i];
                    vw[i] += w[j] * v[indexes[j]][i];
                    xw[i] += w[j] * x[indexes[j]][i];
                }
                uw[i] /= ws;
                vw[i] /= ws;
                xw[i] /= ws;
            }

            R_type s = 0.0;
            R_type t = 0.0;

            // Adapt the covariance matrix and the step size according to Hansen and Ostermeier (2001)
            // and Hansen et al. (2003)
            for (N_type i = 0, i0 = 0; i < n; ++i, i0 += n) {
                pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i]; // ibd. (2001), Eq. (14)
                if (ccov > 0.0) {
                    // BD is not used anymore and can be overwritten
                    valarray<R_type> &Z = BD;

                    for (N_type j = 0, ij = i0; j <= i; ++j, ++ij) {
                        Z[ij] = 0.0;
                        for (N_type k = 0; k < parent_number; ++k) {
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
                decompose(n, C, B, d);

                // Adjust the condition of the covariance matrix and recompute the
                // local step sizes
                if ((t = d[n - 1] / max_covariance_matrix_condition - d[0]) > 0.0) {
                    for (N_type i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                        C[ii] += t;
                        d[i] += t;
                    }
                }
                for (N_type i = 0; i < n; ++i) {
                    d[i] = sqrt(d[i]);
                }
            }

            // Check if the optimization is completed
            for (N_type i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                optimized = (sqr(step_size) * C[ii] <
                             sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
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
    void postopti(const F &f, const Constraint &constraint, N_type n,
                  const R_type x[],
                  const R_type d[],
                  const R_type B[],
                  const R_type C[],
                  R_type s,
                  R_type z[]) {
        using std::abs;
        using std::sqrt;
        using std::valarray;
        using std::thread;

        const R_type zx = f(&x[0], n) + constraint.cost(&x[0], n);

        R_type a = 0.0;
        R_type b = 0.0;
        R_type c = s;

        do {
            // Compute two steps along the line of least variance in opposite directions
            valarray<R_type> p(x, n);
            valarray<R_type> q(x, n);
            for (N_type i = 0, j = 0, ij = j; i < n; ++i, ij += n) {
                p[i] += c * B[ij] * d[j];
                q[i] -= c * B[ij] * d[j];
            }
            R_type zp;
            R_type zq;
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

        for (N_type i = 0, ii = 0; i < n; ++i, ii += n + 1) {
            z[i] = s * sqrt(C[ii]);
        }
    }

}

#endif // ESPECIA_OPTIMIZE_H
