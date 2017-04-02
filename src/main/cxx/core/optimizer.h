/// @file optimizer.h
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
#ifndef ESPECIA_OPTMIZER_H
#define ESPECIA_OPTMIZER_H

#include <cstddef>
#include <functional>

#include "decompose.h"
#include "deviates.h"
#include "mtwister.h"
#include "optimize.h"

using std::valarray;

namespace especia {

    /**
     * No constraint.
     *
     * @tparam T The number type.
     */
    template<class T = R_type>
    class No_Constraint {
    public:
        /**
         * The constructor.
         */
        No_Constraint() {
        }

        /**
         * The destructor.
         */
        ~No_Constraint() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param[in] x The parameter vector.
         * @param[in] n The number of parameters to test.
         * @return always @c false.
         */
        bool is_violated(const T x[], N_type n) const {
            return false;
        }

        /**
         * Computes the cost associated with the constraint.
         *
         * @param[in] x The parameter vector.
         * @param[in] n The number of parameters to take account of.
         * @return always zero.
         */
        T cost(const T x[], N_type n) const {
            return T(0);
        }
    };

    /**
     * No tracing.
     *
     * @tparam T The number type.
     */
    template<class T = R_type>
    class No_Tracing {
    public:

        /**
         * The constructor.
         */
        No_Tracing() {
        }

        /**
         * The destructor.
         */
        ~No_Tracing() {
        }

        /**
         * Tests if tracing is enabled.
         *
         * @param[in] g The generation number.
         * @return always @c false.
         */
        bool is_enabled(L_type g) const {
            return false;
        }

        /**
         * Traces state information.
         *
         * @param[in] g The generation number.
         * @param[in] y The value of the objective function.
         * @param[in] min_step The minimum step size.
         * @param[in] max_step The maximum step size.
         */
        void trace(L_type g, T y, T min_step, T max_step) const {
        }
    };

    /**
     * An optimizer based on the CMA-ES developed by Hansen and Ostermeier (2001).
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
     */
    class Optimizer {
    public:
        /**
         * Builds a new optimizer.
         */
        class Builder {
        public:
            /**
             * Default constructor.
             */
            Builder();

            /**
             * The destructor.
             */
            ~Builder();

            /**
             * Builds a new optimizer.
             *
             * @return the optimizer.
             */
            Optimizer build();

            /**
              * Returns the problem dimension.
              *
              * @return the problem dimension.
              */
            N_type get_problem_dimension() const {
                return n;
            }

            /**
             * Returns the parent number.
             *
             * @return the parent number.
             */
            N_type get_parent_number() const {
                return parent_number;
            }

            /**
             * Returns the population size.
             *
             * @return the population size.
             */
            N_type get_population_size() const {
                return population_size;
            }

            /**
             * Returns the covariance matrix update modulus.
             *
             * @return the covariance matrix update modulus.
             */
            N_type get_covariance_update_modulus() const {
                return update_modulus;
            }

            /**
             * Returns the accuracy goal.
             *
             * @return the accuracy goal.
             */
            R_type get_accuracy_goal() const {
                return accuracy_goal;
            }

            /**
             * Returns the random seed.
             *
             * @return the random seed.
             */
            W_type get_random_seed() const {
                return random_seed;
            }

            /**
             * Returns the stop generation.
             *
             * @return the stop generation.
             */
            L_type get_stop_generation() const {
                return stop_generation;
            }

            /**
             * Returns the recombination weights.
             *
             * @return the recombination weights.
             */
            const std::valarray<R_type> &get_weights() const {
                return weights;
            }

            /**
             * Returns the step size cumulation rate.
             *
             * @return the step size cumulation rate.
             */
            R_type get_step_size_cumulation_rate() const {
                return cs;
            }

            /**
             * Returns the distribution cumulation rate.
             *
             * @return the distribution cumulation rate.
             */
            R_type get_distribution_cumulation_rate() const {
                return cc;
            }

            /**
             * Returns the covariance matrix adaption rate.
             *
             * @return the covariance matrix adaption rate.
             */
            R_type get_covariance_matrix_adaption_rate() const {
                return ccov;
            }

            /**
             * Returns the covariance matrix adaption mixing.
             *
             * @return the covariance matrix adaption mixing.
             */
            R_type get_covariance_matrix_adaption_mixing() const {
                return acov;
            }

            /**
             * Returns the step size damping.
             *
             * @return the step size damping.
             */
            R_type get_step_size_damping() const {
                return step_size_damping;
            }

            /**
             * Configures the problem dimension.
             *
             * @param[in] n The problem dimension.
             * @return this builder.
             */
            Builder &with_problem_dimension(N_type n = 1);

            /**
             * Configures the parent number.
             *
             * @param[in] parent_number The parent number.
             * @return this builder.
             */
            Builder &with_parent_number(N_type parent_number = 4);

            /**
             * Configures the population size.
             *
             * @param[in] population_size The population size.
             * @return this builder.
             */
            Builder &with_population_size(N_type population_size = 8);

            /**
             * Configures the covariance matrix update modulus.
             *
             * @param[in] update_modulus The update modulus.
             * @return this builder.
             */
            Builder &with_covariance_update_modulus(N_type update_modulus = 1);

            /**
             * Configures the accuracy goal.
             *
             * @param[in] accuracy_goal The accuracy goal.
             * @return this builder.
             */
            Builder &with_accuracy_goal(R_type accuracy_goal = 1.0E-04);

            /**
             * Configures the random seed.
             *
             * @param[in] seed The random seed.
             * @return this builder.
             */
            Builder &with_random_seed(W_type seed = 27182);

            /**
             * Configures the stop generation.
             *
             * @param[in] stop_generation The stop generation.
             * @return this builder.
             */
            Builder &with_stop_generation(L_type stop_generation = 1000);

        private:
            /**
             * Returns a pointer to the recombination weights.
             *
             * @return a pointer to the recombination weights.
             */
            const R_type *get_weights_pointer() const {
                return &weights[0];
            }

            /**
             * Calculates strategy parameters like recombination weights, cumulation and adaption rates.
             */
            void set_strategy_parameters();

            /**
             * The problem dimension.
             */
            N_type n = 1;

            /**
             * The parent number.
             */
            N_type parent_number = 4;

            /**
             * The population size.
             */
            N_type population_size = 8;

            /**
             * The covariance matrix update modulus.
             */
            N_type update_modulus = 1;

            /**
             * The accuracy goal.
             */
            R_type accuracy_goal = 1.0E-4;

            /**
              * The random seed.
              */
            W_type random_seed = 27182;

            /**
             * The stop generation.
             */
            L_type stop_generation = 1000;

            /**
             * The recombination weights.
             */
            std::valarray<R_type> weights;

            /**
             * The variance of the recombination weights.
             */
            R_type wv;

            /**
             * The step size cumulation rate.
             */
            R_type cs;

            /**
             * The distribution cumulation rate.
             */
            R_type cc;

            /**
             * The covariance matrix adaption mixing.
             */
            R_type acov;

            /**
             * The covariance matrix adaption rate.
             */
            R_type ccov;

            /**
             * The step size damping.
             */
            R_type step_size_damping;

            friend class Optimizer;
        };

        /**
         * The optimization result.
         */
        class Result {
        public:
            /**
             * The destructor.
             */
            ~Result();

            /**
             * Returns the covariance matrix.
             *
             * @return the covariance matrix.
             */
            const std::valarray<R_type> &get_covariance_matrix() const {
                return C;
            }

            /**
             * Returns the distribution cumulation path.
             *
             * @return the distribution cumulation path.
             */
            const std::valarray<R_type> &get_distribution_cumulation_path() const {
                return pc;
            }

            /**
             * Returns the optimized fitness.
             *
             * @return the optimized fitness.
             */
            R_type get_fitness() const {
                return y;
            }

            /**
             * Returns the final generation number.
             *
             * @return the final generation number.
             */
            L_type get_generation_number() const {
                return g;
            }

            /**
             * Returns the final global step size.
             *
             * @return the final global step size.
             */
            R_type get_global_step_size() const {
                return s;
            }

            /**
             * Returns the final local step sizes.
             *
             * @return the final local step sizes.
             */
            const std::valarray<R_type> &get_local_step_sizes() const {
                return d;
            }

            /**
             * Returns the optimized parameter values.
             *
             * @return the optimized parameter values.
             */
            const std::valarray<R_type> &get_parameter_values() const {
                return x;
            }

            /**
             * Returns the parameter uncertainties.
             *
             * @return the parameter uncertainties.
             */
            const std::valarray<R_type> &get_parameter_uncertainties() const {
                return z;
            }

            /**
             * Returns the final rotation matrix.
             *
             * @return the final rotation matrix.
             */
            const std::valarray<R_type> &get_rotation_matrix() const {
                return B;
            }

            /**
             * Returns the step size cumulation path.
             *
             * @return the step size cumulation path.
             */
            const std::valarray<R_type> &get_step_size_cumulation_path() const {
                return ps;
            }

            /**
             * Returns the optimization status flag.
             *
             * @return the optimization status flag.
             */
            bool is_optimized() const {
                return optimized;
            }

            /**
             * Returns the mutation variance underflow status flag.
             *
             * @return the mutation variance underflow status flag.
             */
            bool is_underflow() const {
                return underflow;
            }

        private:
            /**
             * The constructor.
             *
             * @param[in] n The problem dimension.
             * @param[in] x The initial parameter values.
             * @param[in] d The initial local step sizes.
             * @param[in] s The initial global step size.
             */
            Result(N_type n, const std::valarray<R_type> &x, const std::valarray<R_type> &d, R_type s);

            /**
             * Returns a pointer to the covariance matrix.
             *
             * @return a pointer to the covariance matrix.
             */
            R_type *get_covariance_matrix_pointer() {
                return &C[0];
            }

            /**
             * Returns a pointer to the distribution cumulation path.
             *
             * @return a pointer to the distribution cumulation path.
             */
            R_type *get_distribution_cumulation_path_pointer() {
                return &pc[0];
            }

            /**
             * Returns a reference to the fitness.
             *
             * @return a reference to the fitness.
             */
            R_type &__fitness() {
                return y;
            }

            /**
             * Returns a reference to the generation number.
             *
             * @return a reference to the generation number.
             */
            L_type &__generation_number() {
                return g;
            }

            /**
             * Returns a reference to the global step size.
             *
             * @return a reference to the global step size.
             */
            R_type &__global_step_size() {
                return s;
            }

            /**
             * Returns a pointer to the local step sizes.
             *
             * @return a pointer to the local step sizes.
             */
            R_type *get_local_step_sizes_pointer() {
                return &d[0];
            }

            /**
             * Returns a pointer to the parameter values.
             *
             * @return a pointer to the parameter values.
             */
            R_type *get_parameter_values_pointer() {
                return &x[0];
            }

            /**
             * Returns a pointer to the parameter uncertainties.
             *
             * @return a pointer to the parameter uncertainties.
             */
            R_type *get_parameter_uncertainties_pointer() {
                return &z[0];
            }

            /**
             * Returns a pointer to the rotation matrix.
             *
             * @return a pointer to the rotation matrix.
             */
            R_type *get_rotation_matrix_pointer() {
                return &B[0];
            }

            /**
             * Returns a pointer to the step size cumulation path.
             *
             * @return a pointer to the step size cumulation path.
             */
            R_type *get_step_size_cumulation_path_pointer() {
                return &ps[0];
            }

            /**
             * Returns a reference to the optimization status flag.
             *
             * @return a reference to the optimization status flag.
             */
            bool &__optimized() {
                return optimized;
            }

            /**
             * Returns a reference to the mutation variance underflow status flag.
             *
             * @return a reference to the mutation variance underflow status flag.
             */
            bool &__underflow() {
                return underflow;
            }

            /**
             * The optimized parameter values.
             */
            std::valarray<R_type> x;

            /**
             * The final local step sizes.
             */
            std::valarray<R_type> d;

            /**
             * The final global step size.
             */
            R_type s;

            /**
             * The parameter uncertainties.
             */
            std::valarray<R_type> z;

            /**
             * The optimized fitness.
             */
            R_type y;

            /**
             * The final covariance matrix.
             */
            std::valarray<R_type> C;

            /**
             * The final rotation matrix.
             */
            std::valarray<R_type> B;

            /**
             * The distribution cumulation path.
             */
            std::valarray<R_type> pc;

            /**
             * The step size cumulation path.
             */
            std::valarray<R_type> ps;

            /**
             * The optimization status flag.
             */
            bool optimized;

            /**
             * The mutation variance underflow status flag.
             */
            bool underflow;

            /**
             * The final generation number.
             */
            L_type g;

            friend class Optimizer;
        };

        /**
         * The destructor.
         */
        ~Optimizer();

        /**
         * Maximizes an objective function.
         *
         * @tparam F The function type.
         * @tparam Constraint The constraint type.
         * @tparam Tracing The tracer type.
         *
         * @param[in] f The objective function.
         * @param[in] x The initial parameter values.
         * @param[in] d The initial local step sizes.
         * @param[in] s The initial global step size.
         * @param[in] constraint The constraint.
         * @param[in] tracer The tracer.
         *
         * @return the maximization result.
         */
        template<class F, class Constraint, class Tracing>
        Result maximize(const F &f,
                        const std::valarray<R_type> &x,
                        const std::valarray<R_type> &d,
                        const R_type s,
                        const Constraint &constraint = No_Constraint<>(),
                        const Tracing &tracer = No_Tracing<>()) const {
            return optimize(f, x, d, s, constraint, tracer, std::greater<R_type>());
        }

        /**
         * Minimizes an objective function.
         *
         * @tparam F The function type.
         * @tparam Constraint The constraint type.
         * @tparam Tracing The tracer type.
         *
         * @param[in] f The objective function.
         * @param[in] x The initial parameter values.
         * @param[in] d The initial local step sizes.
         * @param[in] s The initial global step size.
         * @param[in] constraint The constraint.
         * @param[in] tracer The tracer.
         *
         * @return the minimization result.
         */
        template<class F, class Constraint, class Tracing>
        Result minimize(const F &f,
                        const std::valarray<R_type> &x,
                        const std::valarray<R_type> &d,
                        const R_type s,
                        const Constraint &constraint = No_Constraint<>(),
                        const Tracing &tracer = No_Tracing<>()) const {
            return optimize(f, x, d, s, constraint, tracer, std::less<R_type>());
        }

    private:
        /**
         * Creates a new instance of this class with the build configuration supplied as argument.
         *
         * @param[in] builder The build configuration.
         */
        Optimizer(const Builder &builder);

        /**
         * Optimizes an objective function.
         *
         * @tparam F The function type.
         * @tparam Constraint The constraint type.
         * @tparam Tracing The tracer type.
         * @tparam Compare The fitness comparator type.
         *
         * @param[in] f The objective function.
         * @param[in] x The initial parameter values.
         * @param[in] d The initial local step sizes.
         * @param[in] s The initial global step size.
         * @param[in] constraint The constraint.
         * @param[in] tracer The tracer.
         * @param[in] compare The fitness comparator.
         *
         * @return the optimization result.
         */
        template<class F, class Constraint, class Tracing, class Compare>
        Result optimize(const F &f,
                        const std::valarray<R_type> &x,
                        const std::valarray<R_type> &d,
                        const R_type s,
                        const Constraint &constraint,
                        const Tracing &tracer,
                        const Compare &compare) const {
            using especia::optimize;
            using especia::postopti;

            const N_type n = config.get_problem_dimension();

            Result result(n, x, d, s);

            optimize(f, constraint, n,
                     config.get_parent_number(),
                     config.get_population_size(),
                     config.get_weights_pointer(),
                     config.get_step_size_damping(),
                     config.get_step_size_cumulation_rate(),
                     config.get_distribution_cumulation_rate(),
                     config.get_covariance_matrix_adaption_rate(),
                     config.get_covariance_matrix_adaption_mixing(),
                     config.get_covariance_update_modulus(),
                     config.get_accuracy_goal(),
                     config.get_stop_generation(),
                     result.__generation_number(),
                     result.get_parameter_values_pointer(),
                     result.__global_step_size(),
                     result.get_local_step_sizes_pointer(),
                     result.get_rotation_matrix_pointer(),
                     result.get_covariance_matrix_pointer(),
                     result.get_step_size_cumulation_path_pointer(),
                     result.get_distribution_cumulation_path_pointer(),
                     result.__fitness(),
                     result.__optimized(),
                     result.__underflow(),
                     deviate, decompose, compare, tracer
            );

            if (result.__optimized()) {
                postopti(f, constraint, n,
                         result.get_parameter_values_pointer(),
                         result.get_local_step_sizes_pointer(),
                         result.get_rotation_matrix_pointer(),
                         result.get_covariance_matrix_pointer(),
                         result.get_global_step_size(),
                         result.get_parameter_uncertainties_pointer()
                );
            }

            return result;
        }

        /**
         * The build configuration.
         */
        const Builder config;

        /**
         * The eigenvalue decomposition strategy.
         */
        const Decompose decompose;

        /**
         * The random number generator.
         */
        const Normal_Deviate<Mt19937> deviate;
    };

}

#endif // ESPECIA_OPTMIZER_H
