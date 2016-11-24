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
#ifndef RQ_OPTIMIZE_H
#define RQ_OPTIMIZE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <valarray>

namespace RQ {
    // CMA-ES function template for parametric nonlinear function optimization,
    // based on Hansen and Ostermeier (2001).
    template<class objective_function, class normal_deviate, class sym_eig_decomp, class comparation>
    void optimize(objective_function &f, double x[], size_t n,
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // CMA-ES function template for parametric nonlinear function optimization,
    // based on Hansen and Ostermeier (2001).
    template<class objectp, class functionp, class normal_deviate, class sym_eig_decomp,
            class comparation>
    void optimize(objectp obj, functionp f, double x[], size_t n,
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // CMA-ES function template for constrained (hard box) nonlinear function
    // optimization, based on Hansen and Ostermeier (2001).
    template<class objective_function, class normal_deviate, class sym_eig_decomp, class comparation>
    void optimize(objective_function &f, double x[], size_t n,
                  const double inf[], const double sup[],
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // CMA-ES function template for constrained (hard box) nonlinear function
    // optimization, based on Hansen and Ostermeier (2001).
    template<class objectp, class functionp, class normal_deviate, class sym_eig_decomp,
            class comparation>
    void optimize(objectp obj, functionp f, double x[], size_t n,
                  const double inf[], const double sup[],
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // CMA-ES function template for constrained nonlinear function optimization,
    // based on Hansen and Ostermeier (2001).
    template<class objective_function, class constraint, class normal_deviate,
            class sym_eig_decomp, class comparation>
    void optimize(objective_function &f, double x[], size_t n, constraint &reject,
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // CMA-ES function template for constrained nonlinear function optimization,
    // based on Hansen and Ostermeier (2001).
    template<class objectp, class functionp, class constraint, class normal_deviate,
            class sym_eig_decomp, class comparation>
    void optimize(objectp obj, functionp f, double x[], size_t n, constraint &reject,
                  size_t parent_number,
                  size_t population_size,
            // twice the parent number, at least
                  const double weight[],
                  double &step_size,
                  double step_size_damping,
                  double step_size_cumulation,
                  double distribution_cumulation,
                  double covariance_matrix_adaption_rate,
                  double covariance_matrix_adaption_mixing,
                  double accuracy_goal,
                  unsigned long stop_generation,
                  double diagonal_matrix[],
            // diagonal matrix (packed)
                  double rotation_matrix[],
            // orthogonal matrix (row-major)
                  double covariance_matrix[],
            // symmetric matrix (row-major, lower triangular = column-major, upper triangular)
                  double step_size_evolution_path[],
                  double distribution_evolution_path[],
                  unsigned long &generation_number,
                  bool &is_optimized,
                  bool &mutation_variance_underflow,
                  double &optimum,
                  normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned update_modulus = 1);

    // Function template to calculate uncertainty and covariance
    template<class objective_function>
    void scale_step_size(objective_function &f, const double x[], size_t n,
                         double &step_size,
                         const double diagonal_matrix[],
            // diagonal matrix (packed)
                         const double rotation_matrix[]
            // orthogonal matrix (row-major)
    );

    // Function template to calculate uncertainty and covariance
    template<class objectp, class functionp>
    void scale_step_size(objectp obj, functionp f, const double x[], size_t n,
                         double &step_size,
                         const double diagonal_matrix[],
            // diagonal matrix (packed)
                         const double rotation_matrix[]
            // orthogonal matrix (row-major)
    );

    template<class number>
    number sqr(number x);

    template<class number, class comparation>
    class indirect_comparation;
}

template<class objective_function, class normal_deviate, class sym_eig_decomp, class comparation>
void
RQ::optimize(objective_function &f, double xw[], size_t n,
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
    const double csu = sqrt(cs * (2.0 - cs));
    const double ccu = sqrt(cc * (2.0 - cc));
    const double ws = accumulate(w, w + parent_number, 0.0);
    const double cw = ws / sqrt(inner_product(w, w + parent_number, w, 0.0));

    while (g < stop_generation) {
        valarray<double> uw(n);
        valarray<double> zw(n);
        valarray<valarray<double> > u(uw, population_size);
        valarray<valarray<double> > x = u;
        valarray<valarray<double> > z = u;

        valarray<double> fitness(population_size);
        valarray<size_t> index(population_size);

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness
        for (size_t k = 0; k < population_size; ++k) {
            for (size_t i = 0; i < n; ++i)
                z[k][i] = ndev();
            for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
                u[k][i] = inner_product(&z[k][0], &z[k][n], &BD[i0], 0.0);
                x[k][i] = xw[i] + step_size * u[k][i];
                // Hansen and Ostermeier (2001), Eq. (13)
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
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
            break;

        // Recombine the best individuals
        for (size_t i = 0; i < n; ++i) {
            uw[i] = xw[i] = zw[i] = 0.0;
            for (size_t j = 0; j < parent_number; ++j) {
                uw[i] += w[j] * u[index[j]][i];
                xw[i] += w[j] * x[index[j]][i];
                zw[i] += w[j] * z[index[j]][i];
            }
            uw[i] /= ws;
            xw[i] /= ws;
            zw[i] /= ws;
        }

        double s = 0.0;
        double t = 0.0;

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * inner_product(&zw[0], &zw[n], &B[i0], 0.0);
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = f(xw, n);
}

template<class objectp, class functionp, class normal_deviate, class sym_eig_decomp, class comparation>
void
RQ::optimize(objectp obj, functionp f, double xw[], size_t n,
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
    const double csu = sqrt(cs * (2.0 - cs));
    const double ccu = sqrt(cc * (2.0 - cc));
    const double ws = accumulate(w, w + parent_number, 0.0);
    const double cw = ws / sqrt(inner_product(w, w + parent_number, w, 0.0));

    while (g < stop_generation) {
        valarray<double> uw(n);
        valarray<double> zw(n);
        valarray<valarray<double> > u(uw, population_size);
        valarray<valarray<double> > x = u;
        valarray<valarray<double> > z = u;

        valarray<double> fitness(population_size);
        valarray<size_t> index(population_size);

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness
        for (size_t k = 0; k < population_size; ++k) {
            for (size_t i = 0; i < n; ++i)
                z[k][i] = ndev();
            for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
                u[k][i] = inner_product(&z[k][0], &z[k][n], &BD[i0], 0.0);
                x[k][i] = xw[i] + step_size * u[k][i];
                // Hansen and Ostermeier (2001), Eq. (13)
            }
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t k = 0; k < population_size; ++k) {
            fitness[k] = (obj->*f)(&x[k][0], n);
            index[k] = k;
        }
        partial_sort(&index[0], &index[parent_number], &index[population_size],
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
            break;

        // Recombine the best individuals
        for (size_t i = 0; i < n; ++i) {
            uw[i] = xw[i] = zw[i] = 0.0;
            for (size_t j = 0; j < parent_number; ++j) {
                uw[i] += w[j] * u[index[j]][i];
                xw[i] += w[j] * x[index[j]][i];
                zw[i] += w[j] * z[index[j]][i];
            }
            uw[i] /= ws;
            xw[i] /= ws;
            zw[i] /= ws;
        }

        double s = 0.0;
        double t = 0.0;

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * inner_product(&zw[0], &zw[n], &B[i0], 0.0);
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = (obj->*f)(xw, n);
}

template<class objective_function, class normal_deviate, class sym_eig_decomp, class comparation>
void
RQ::optimize(objective_function &f, double xw[], size_t n,
             const double inf[], const double sup[],
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
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

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness        
        for (size_t k = 0; k < population_size; ++k) {
            uw = vw = 0.0;
            // buffer, initialization is essential
            for (size_t j = 0, reject; j < n; ++j) {
                do {
                    const double z = ndev();

                    for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                        u[k][i] = uw[i] + z * BD[ij];
                        v[k][i] = vw[i] + z * B[ij];
                        x[k][i] = xw[i] + u[k][i] * step_size;
                        // Hansen and Ostermeier (2001), Eq. (13)

                        reject = (x[k][i] < inf[i] or x[k][i] > sup[i]);
                        if (reject)
                            break;
                    }
                } while (reject);

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
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
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

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i];
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = f(xw, n);
}

template<class objectp, class functionp, class normal_deviate, class sym_eig_decomp, class comparation>
void
RQ::optimize(objectp obj, functionp f, double xw[], size_t n,
             const double inf[], const double sup[],
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
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

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness        
        for (size_t k = 0; k < population_size; ++k) {
            uw = vw = 0.0;
            // buffer, initialization is essential
            for (size_t j = 0, reject; j < n; ++j) {
                do {
                    const double z = ndev();

                    for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                        u[k][i] = uw[i] + z * BD[ij];
                        v[k][i] = vw[i] + z * B[ij];
                        x[k][i] = xw[i] + u[k][i] * step_size;
                        // Hansen and Ostermeier (2001), Eq. (13)

                        reject = (x[k][i] < inf[i] or x[k][i] > sup[i]);
                        if (reject)
                            break;
                    }
                } while (reject);

                uw = u[k];
                vw = v[k];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t k = 0; k < population_size; ++k) {
            fitness[k] = (obj->*f)(&x[k][0], n);
            index[k] = k;
        }
        partial_sort(&index[0], &index[parent_number], &index[population_size],
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
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

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i];
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = (obj->*f)(xw, n);
}

template<class objective_function, class constraint, class normal_deviate,
        class sym_eig_decomp, class comparation>
void
RQ::optimize(objective_function &f, double xw[], size_t n, constraint &reject,
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
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

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness        
        for (size_t k = 0; k < population_size; ++k) {
            uw = vw = 0.0;
            // buffer, initialization is essential
            for (size_t j = 0; j < n; ++j) {
                do {
                    const double z = ndev();

                    for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                        u[k][i] = uw[i] + z * BD[ij];
                        v[k][i] = vw[i] + z * B[ij];
                        x[k][i] = xw[i] + u[k][i] * step_size;
                        // Hansen and Ostermeier (2001), Eq. (13)
                    }
                } while (reject);

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
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
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

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i];
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = f(xw, n);
}

template<class objectp, class functionp, class constraint, class normal_deviate,
        class sym_eig_decomp, class comparation>
void
RQ::optimize(objectp obj, functionp f, double xw[], size_t n, constraint &reject,
             size_t parent_number,
             size_t population_size,
             const double w[],
             double &step_size,
             double step_size_damping,
             double cs,
             double cc,
             double ccov,
             double acov,
             double accuracy_goal,
             unsigned long stop_generation,
             double d[],
             double B[],
             double C[],
             double ps[],
             double pc[],
             unsigned long &g,
             bool &is_opt,
             bool &is_ufl,
             double &optimum,
             normal_deviate &ndev, sym_eig_decomp &evd, comparation comp, unsigned um) {
    using std::accumulate;
    using std::exp;
    using std::inner_product;
    using std::numeric_limits;
    using std::partial_sort;
    using std::sqrt;
    using std::valarray;

    const double expected_length = (n - 0.25 + 1.0 / (21 * n)) / sqrt(double(n));
    const double max_covariance_matrix_condition = 0.01 /
                                                   numeric_limits<double>::epsilon();
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

        // Buffer BD
        valarray<double> BD(B, n * n);
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0, ij = j; i < n; ++i, ij += n)
                BD[ij] *= d[j];

        // Generate a new population of object parameter vectors,
        // sorted indirectly by fitness        
        for (size_t k = 0; k < population_size; ++k) {
            uw = vw = 0.0;
            // buffer, initialization is essential
            for (size_t j = 0; j < n; ++j) {
                do {
                    const double z = ndev();

                    for (size_t i = 0, ij = j; i < n; ++i, ij += n) {
                        u[k][i] = uw[i] + z * BD[ij];
                        v[k][i] = vw[i] + z * B[ij];
                        x[k][i] = xw[i] + u[k][i] * step_size;
                        // Hansen and Ostermeier (2001), Eq. (13)
                    }
                } while (reject);

                uw = u[k];
                vw = v[k];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t k = 0; k < population_size; ++k) {
            fitness[k] = (obj->*f)(&x[k][0], n);
            index[k] = k;
        }
        partial_sort(&index[0], &index[parent_number], &index[population_size],
                     indirect_comparation<double, comparation>(fitness, comp));
        ++g;

        // Check the mutation variance
        is_ufl = (fitness[index[0]] == fitness[index[parent_number]]);
        if (!is_ufl)
            for (size_t i = 0, ij = g % n; i < n; ++i, ij += n) {
                is_ufl = (xw[i] == xw[i] + 0.2 * step_size * BD[ij]);
                if (!is_ufl)
                    break;
            }
        if (is_ufl)
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

        // Adapt the covariance matrix and the step size according to
        // Hansen and Ostermeier (2001)
        for (size_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
            pc[i] = (1.0 - cc) * pc[i] + (ccu * cw) * uw[i];
            // ibd., Eq. (14)
            if (ccov > 0.0) {
                valarray<double> &Z = BD;
                // BD is not used anymore, can be overwritten
                for (size_t j = 0, ij = i0; j <= i; ++j, ++ij) {
                    Z[ij] = 0.0;
                    for (size_t k = 0; k < parent_number; ++k)
                        Z[ij] += w[k] * (u[index[k]][i] * u[index[k]][j]);
                    C[ij] = (1.0 - ccov) * C[ij] + ccov * (acov * (pc[i] * pc[j]) +
                                                           (1.0 - acov) * Z[ij] / ws);
                    // Hansen et al. (2003), Eq. (11)
                }
            }
            ps[i] = (1.0 - cs) * ps[i] + (csu * cw) * vw[i];
            // ibd., Eq. (16)
            s += ps[i] * ps[i];
        }
        step_size *= exp((cs / step_size_damping) * (sqrt(s) / expected_length - 1.0));
        // ibd., Eq. (17)

        if (ccov > 0.0 and g % um == 0) {
            // Decompose the covariance matrix and sort its eigenvalues in ascending
            // order, along with eigenvectors
            evd(C, B, d, n);

            // Adjust the condition of the covariance matrix and compute the
            // square roots of its eigenvalues
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
            is_opt = (sqr(step_size) * C[ii] < sqr(accuracy_goal * xw[i]) + 1.0 / max_covariance_matrix_condition);
            if (!is_opt)
                break;
        }
        if (is_opt)
            break;
    }

    optimum = (obj->*f)(xw, n);
}

template<class objective_function>
void
RQ::scale_step_size(objective_function &f, const double x[], size_t n,
                    double &s,
                    const double d[],
                    const double B[]) {
    using std::abs;
    using std::valarray;

    const double a = 100.0 * s;

    valarray<double> p(x, n);
    valarray<double> q(x, n);
    for (size_t i = 0, j = n - 1, ij = j; i < n; ++i, ij += n) {
        p[i] += a * B[ij] * d[n - 1];
        q[i] += a * B[ij] * d[n - 1];
    }

    const double zx = f(&x[0], n);
    const double zp = f(&p[0], n);
    const double zq = f(&q[0], n);
    s = a / sqrt(abs(2.0 * (zp - zx) - (zp - zq)));
    // compute the covariance along the major principal axis by means of a parabola
}

template<class objectp, class functionp>
void
RQ::scale_step_size(objectp obj, functionp f, const double x[], size_t n,
                    double &s,
                    const double d[],
                    const double B[]) {
    using std::abs;
    using std::valarray;

    const double a = 100.0 * s;

    valarray<double> p(x, n);
    valarray<double> q(x, n);
    for (size_t i = 0, j = n - 1, ij = j; i < n; ++i, ij += n) {
        p[i] += a * B[ij] * d[j];
        q[i] += a * B[ij] * d[j];
    }

    const double zx = (obj->*f)(&x[0], n);
    const double zp = (obj->*f)(&p[0], n);
    const double zq = (obj->*f)(&q[0], n);
    s = a / sqrt(abs(2.0 * (zp - zx) - (zp - zq)));
    // compute the covariance along the major principal axis by means of a parabola
}

template<class number>
inline
number
RQ::sqr(number x) {
    return (x == number(0)) ? number(0) : x * x;
}

template<class number, class comparation>
class RQ::indirect_comparation {
public:
    indirect_comparation(const std::valarray<number> &fitness, const comparation &comp);

    ~indirect_comparation();

    bool operator()(size_t i, size_t j);
    // comparation operator

private:
    const std::valarray<number> &fitness;
    const comparation &comp;
};

template<class number, class comparation>
RQ::indirect_comparation<number, comparation>::indirect_comparation(const std::valarray<number> &f,
                                                                    const comparation &c) : fitness(f), comp(c) {
}

template<class number, class comparation>
RQ::indirect_comparation<number, comparation>::~indirect_comparation() {
}

template<class number, class comparation>
inline
bool
RQ::indirect_comparation<number, comparation>::operator()(size_t i, size_t j) {
    return comp(fitness[i], fitness[j]);
}

#endif // RQ_OPTIMIZE_H

// References
//
// N. Hansen, S. D. M"uller, P. Koumoutsakos (2003)
//   Reducing the Increasing the Time Complexity of the Derandomized Evolution
//     Strategy with Covariance Matrix Adaption (CMA-ES)
//   Evolutionary Computation, 11, 1, ISSN 1063-6560
//
// N. Hansen, A. Ostermeier (2001)
//   Completely Derandomized Self-Adaption in Evolution Strategies
//   Evolutionary Computation, 9, 159, ISSN 1063-6560
