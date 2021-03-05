/// @file optimize.h
/// CMA-ES function templates for nonlinear function optimization.
/// Copyright (c) 2021 Ralf Quast
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
        ~Index_Compare() = default;

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
     * Evolution strategy with covariance matrix adaption (CMA-ES) for nonlinear function optimization.
     * Based on Hansen (2014, http://cma.gforge.inria.fr/purecmaes.m).
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
     * @param[in] ccov The rank-1 covariance matrix adaption rate.
     * @param[in] acov The rank-µ covariance matrix adaption rate.
     * @param[in] update_modulus The covariance matrix update modulus.
     * @param[in] accuracy_goal The accuracy goal.
     * @param[in] stop_generation The stop generation.
     * @param[in,out] g The generation number.
     * @param[in,out] xw The parameter values.
     * @param[in,out] step_size The global step size.
     * @param[in,out] d The local step sizes.
     * @param[in,out] B The rotation matrix (in column-major layout).
     * @param[in,out] C The covariance matrix (upper triangular part only, in column-major layout).
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
        using std::numeric_limits;
        using std::partial_sort;
        using std::sqrt;
        using std::thread;
        using std::valarray;
        using std::vector;

        const real expected_norm = (n - 0.25 + 1.0 / (21 * n)) / sqrt(real(n));
        const real max_covariance_matrix_condition = 0.01 / numeric_limits<real>::epsilon();
        const real csu = sqrt(cs * (2.0 - cs));
        const real ccu = sqrt(cc * (2.0 - cc));
        const real ws = accumulate(w, w + parent_number, 0.0);
        const real cw = ws / norm(parent_number, w);

        valarray<real> uw(n);
        valarray<real> vw(n);
        valarray<valarray<real>> u(uw, population_size);
        valarray<valarray<real>> v = u;
        valarray<valarray<real>> x = u;

        valarray<real> y(population_size);
        valarray<natural> indexes(population_size);

        while (g < stop_generation) {
            // Generate a new population of object parameter vectors,
            // sorted indirectly by fitness
            for (natural k = 0; k < population_size; ++k) {
                uw = 0.0;
                vw = 0.0;
                for (natural j = 0, nj = 0; j < n; ++j, nj += n) {
                    do {
                        const real z = deviate();

                        for (natural i = 0, ij = nj; i < n; ++i, ++ij) {
                            u[k][i] = uw[i] + z * (B[ij] * d[j]);
                            v[k][i] = vw[i] + z * B[ij];
                            x[k][i] = xw[i] + u[k][i] * step_size; // Hansen & Ostermeier (2001, Eq. 13)
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
            if (underflow) {
                break;
            }

            // Recombine the best individuals
            for (natural i = 0; i < n; ++i) {
                uw[i] = vw[i] = xw[i] = 0.0;
                for (natural k = 0; k < parent_number; ++k) {
                    uw[i] += w[k] * u[indexes[k]][i];
                    vw[i] += w[k] * v[indexes[k]][i];
                    xw[i] += w[k] * x[indexes[k]][i];
                }
                uw[i] /= ws;
                vw[i] /= ws;
                xw[i] /= ws;
            }

            // Adapt the covariance matrix and the step size according to Hansen & Ostermeier (2001)
            // and Hansen (2014)
            if (acov > 0.0 or ccov > 0.0) {
                for (natural j = 0, nj = 0; j < n; ++j, nj += n) {
                    pc[j] = (1.0 - cc) * pc[j] + (ccu * cw) * uw[j]; // Hansen & Ostermeier (2001, Eq. 14)
                    for (natural i = 0, ij = nj; i <= j; ++i, ++ij) {
                        real z = 0.0;
                        for (natural k = 0; k < parent_number; ++k) {
                            z += w[k] * (u[indexes[k]][i] * u[indexes[k]][j]);
                        }
                        // Hansen (2014, http://www.lri.fr/~hansen/purecmaes.m)
                        C[ij] = (C[ij] + acov * (pc[i] * pc[j] - C[ij])) + ccov * (z / ws - C[ij]);
                    }
                }
                if (g % update_modulus == 0) {
                    decompose(C, B, d);

                    const real t = d[n - 1] / max_covariance_matrix_condition - d[0];
                    if (t > 0.0) {
                        for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                            C[ii] += t;
                            d[i] += t;
                        }
                    }
                    for (natural i = 0; i < n; ++i) {
                        d[i] = sqrt(d[i]);
                    }
                }
            }
            for (natural i = 0; i < n; ++i) {
                ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i]; // Hansen & Ostermeier (2001, Eq. 16)
            }
            // Hansen & Ostermeier (2001, Eq. 17)
            step_size *= exp((cs / step_size_damping) * (norm(n, ps) / expected_norm - 1.0));

            // Check if the optimization is completed
            for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
                optimized = (sq(step_size) * C[ii] < sq(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
                if (!optimized) {
                    break;
                }
            }
            if (optimized or tracer.is_tracing(g)) {
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
     * Computes the standard variance along ellipsoid principal axes from the curvature of a
     * parabola through three points around the minimum. The global step size is rescaled to
     * approximate the standard covariance matrix.
     *
     * @tparam F The function type.
     * @tparam Constraint The constraint type.
     *
     * @param[in] f The objective function.
     * @param[in] constraint The constraint on parameter values.
     * @param[in] n The number of parameter values.
     * @param[in] x The parameter values.
     * @param[in] d The local step sizes
     * @param[in] B The rotation matrix (in column-major layout).
     * @param[in] C The covariance matrix (upper triangular part only, in column-major layout).
     * @param[in] s The global step size.
     * @param[out] z The parameter uncertainties.
     */
    template<class F, class Constraint>
    void postopti(const F &f, const Constraint &constraint, natural n,
                  const real x[],
                  const real d[],
                  const real B[],
                  const real C[],
                  const real s,
                  real z[]) {
        using std::abs;
        using std::exp;
        using std::log;
        using std::sqrt;
        using std::valarray;
        using std::thread;

        const real zx = f(&x[0], n) + constraint.cost(&x[0], n);
        // The rescaled global step sizes
        valarray<real> g(s, n);
        
        for (natural j = 0; j < n; ++j) {
            real a = 0.0;
            real b = 0.0;
            real c = g[j];
            
            do {
                // Compute two steps along a principal axis in opposite directions
                valarray<real> p(x, n);
                valarray<real> q(x, n);
                for (natural i = 0, ij = j * n; i < n; ++i, ++ij) {
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
                g[j] = c / sqrt(abs((zp + zq) - (zx + zx)));

                // Make a smaller or larger computation step in the next iteration
                if (abs(0.5 * (zp + zq) - zx) < 0.5) {
                    a = c;
                    c = c * 1.618;
                } else {
                    b = c;
                    c = c * 0.618;
                }
            } while (a == 0.0 or b == 0.0); // the computation step is too small or too large
        }
        // Take the geometric mean to rescale the covariance matrix
        const real h = exp(g.apply(log).sum() / real(n));
        
        for (natural i = 0, ii = 0; i < n; ++i, ii += n + 1) {
            z[i] = h * sqrt(C[ii]);
        }
    }

}

#endif // ESPECIA_OPTIMIZE_H
