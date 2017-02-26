// CMA-ES classes for nonlinear function optimization
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

#ifndef ESPECIA_OPTMIZER_H
#define ESPECIA_OPTMIZER_H

#include <cstddef>
#include <iomanip>
#include <iostream>

#include "optimize.h"

namespace especia {

    /**
     * A bounded constraint.
     *
     * @tparam Number The number type.
     */
    template<class Number>
    class Bounds {
    public:
        /**
         * Constructs a new strict-bound prior constraint.
         *
         * @param lower_bounds[in] The lower bounds.
         * @param upper_bounds[in] The upper bounds.
         * @param n The number of bounds.
         */
        Bounds(const Number lower_bounds[], const Number upper_bounds[], size_t n)
                : a(lower_bounds, n), b(upper_bounds, n) {
        }

        /**
         * Destructor.
         */
        ~Bounds() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param x[in] The parameter vector.
         * @param n[in] The number of parameters to test.
         * @return @c true, if the parameter vector violates the constraint.
         */
        bool is_violated(const Number x[], size_t n) const {
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
         * @param x[in] The parameter vector.
         * @param n[in] The number of parameters to take account of.
         * @return always zero.
         */
        Number cost(const Number x[], size_t n) const {
            return Number(0);
        }

    private:
        const std::valarray<Number> a;
        const std::valarray<Number> b;
    };

    /**
     * No constraint.
     *
     * @tparam Number The number type.
     */
    template<class Number>
    class Unconstrained {
    public:
        /**
         * Constructor.
         */
        Unconstrained() {
        }

        /**
         * Destructor.
         */
        ~Unconstrained() {
        }

        /**
         * Tests if a given parameter vector violates the constraint.
         *
         * @param x[in] The parameter vector.
         * @param n[in] The number of parameters to test.
         * @return always @c false.
         */
        bool is_violated(const Number *x, size_t n) const {
            return false;
        }

        /**
         * Computes the cost associated with the constraint.
         *
         * @param x[in] The parameter vector.
         * @param n[in] The number of parameters to take account of.
         * @return always zero.
         */
        Number cost(const Number x[], size_t n) const {
            return Number(0);
        }
    };

    /**
     * Traces state information to an output stream.
     *
     * @tparam Number The number type.
     */
    template<class Number>
    class Output_Stream_Tracer {
    public:
        /**
         * Constructor.
         *
         * @param output_stream[in] The output stream.
         * @param modulus[in] The trace modulus.
         * @param precision[in] The precision of numeric output.
         * @param width[in] The width of the numeric output fields.
         */
        Output_Stream_Tracer(std::ostream &output_stream, unsigned int modulus, unsigned int precision = 4,
                             unsigned int width = 12)
                : os(output_stream), m(modulus), p(precision), w(width) {
        }

        /**
         * Destructor.
         */
        ~Output_Stream_Tracer() {
        }

        /**
         * Tests if tracing is enabled.
         *
         * @param g[in] The generation number.
         * @return @true if tracing is enabled, otherwise @c false.
         */
        bool is_enabled(unsigned long g) const {
            return m > 0 and g % m == 0;
        }

        /**
         * Traces state information to an output stream..
         *
         * @param g[in] The generation number.
         * @param y[in] The value of the objective function.
         * @param min_step[in] The minimum step size.
         * @param max_step[in] The maximum step size.
         */
        void trace(unsigned long g, Number y, Number min_step, Number max_step) const {
            using std::endl;
            using std::ios_base;
            using std::setw;

            const ios_base::fmtflags fmt = os.flags();

            os.setf(ios_base::fmtflags());
            os.setf(ios_base::scientific, ios_base::floatfield);
            os.setf(ios_base::right, ios_base::adjustfield);
            os.precision(p);

            os << setw(8) << g;
            os << setw(w) << y;
            os << setw(w) << min_step;
            os << setw(w) << max_step;
            os << endl;

            os.flags(fmt);
        }

    private:
        std::ostream &os;
        const unsigned int m;
        const unsigned int p;
        const unsigned int w;
    };


    /**
     * No tracing.
     *
     * @tparam Number The number type.
     */
    template<class Number>
    class No_Tracing {
    public:

        /**
         * Constructor.
         */
        No_Tracing() {
        }

        /**
         * Destructor.
         */
        ~No_Tracing() {
        }

        /**
         * Tests if tracing is enabled.
         *
         * @param g[in] The generation number.
         * @return always @c false.
         */
        bool is_enabled(unsigned long g) const {
            return false;
        }

        /**
         * Traces state information.
         *
         * @param g[in] The generation number.
         * @param y[in] The value of the objective function.
         * @param min_step[in] The minimum step size.
         * @param max_step[in] The maximum step size.
         */
        void trace(unsigned long g, Number y, Number min_step, Number max_step) const {
        }
    };

    /**
     * An optimizer based on the CMA-ES.
     *
     * @todo - this is work in progress
     */
    class Optimizer {
    public:

        /**
         * Builds a new optimizer.
         */
        class Builder {
        public:
            /**
             * Creates a new instance of this class.
             */
            Builder();

            /**
             * Destructor.
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
            size_t get_problem_dimension() const {
                return n;
            }

            /**
             * Returns the parent number.
             *
             * @return the parent number.
             */
            unsigned int get_parent_number() const {
                return parent_number;
            }

            /**
             * Returns the population size.
             *
             * @return the population size.
             */
            unsigned int get_population_size() const {
                return population_size;
            }

            /**
             * Returns the covariance matrix update modulus.
             *
             * @return the covariance matrix update modulus.
             */
            unsigned int get_update_modulus() const {
                return update_modulus;
            }

            /**
             * Returns the accuracy goal.
             *
             * @return the accuracy goal.
             */
            double get_accuracy_goal() const {
                return accuracy_goal;
            }

            /**
             * Returns the stop generation.
             *
             * @return the stop generation.
             */
            unsigned long get_stop_generation() const {
                return stop_generation;
            }

            /**
             * Returns the recombination weights.
             *
             * @return the recombination weights.
             */
            std::valarray<double> get_weights() const {
                return weights;
            }

            /**
             * Returns the step size cumulation rate.
             *
             * @return the step size cumulation rate.
             */
            double get_step_size_cumulation_rate() const {
                return cs;
            }

            /**
             * Returns the distribution cumulation rate.
             *
             * @return the distribution cumulation rate.
             */
            double get_distribution_cumulation_rate() const {
                return cc;
            }

            /**
             * Returns the covariance matrix adaption rate.
             *
             * @return the covariance matrix adaption rate.
             */
            double get_covariance_matrix_adaption_rate() const {
                return ccov;
            }

            /**
             * Returns the covariance matrix adaption mixing.
             *
             * @return the covariance matrix adaption mixing.
             */
            double get_covariance_matrix_adaption_mixing() const {
                return acov;
            }

            /**
             * Returns the step size damping.
             *
             * @return the step size damping.
             */
            double get_step_size_damping() const {
                return step_size_damping;
            }

            /**
             * Configures the problem dimension.
             *
             * @param n[in] The problem dimension.
             * @return this builder.
             */
            Builder &with_problem_dimension(unsigned n = 1);

            /**
             * Configures the parent number.
             *
             * @param parent_number[in] The parent number.
             * @return this builder.
             */
            Builder &with_parent_number(unsigned parent_number = 4);

            /**
             * Configures the population size.
             *
             * @param population_size[in] The population size.
             * @return this builder.
             */
            Builder &with_population_size(unsigned population_size = 8);

            /**
             * Configures the covariance matrix update modulus.
             *
             * @param update_modulus[in] The update modulus.
             * @return this builder.
             */
            Builder &with_update_modulus(unsigned update_modulus = 1);

            /**
             * Configures the accuracy goal.
             *
             * @param accuracy_goal[in] The accuracy goal.
             * @return this builder.
             */
            Builder &with_accuracy_goal(double accuracy_goal = 1.0E-04);

            /**
             * Configures the stop generation.
             *
             * @param stop_generation[in] The stop generation.
             * @return this builder.
             */
            Builder &with_stop_generation(unsigned long stop_generation = 1000);

        private:
            /**
             * Calculates strategy parameters like recombination weights, cumulation and adaption rates.
             */
            void set_strategy_parameters();

            /**
             * The problem dimension.
             */
            size_t n;

            /**
             * The parent number.
             */
            unsigned int parent_number;

            /**
             * The population size.
             */
            unsigned int population_size;

            /**
             * The covariance matrix update modulus.
             */
            unsigned int update_modulus;

            /**
             * The accuracy goal.
             */
            double accuracy_goal;

            /**
             * The stop generation.
             */
            unsigned long stop_generation;

            /**
             * The recombination weights.
             */
            std::valarray<double> weights;

            /**
             * The variance of the recombination weights.
             */
            double wv;

            /**
             * The step size cumulation rate.
             */
            double cs;

            /**
             * The distribution cumulation rate.
             */
            double cc;

            /**
             * The covariance matrix adaption mixing.
             */
            double acov;

            /**
             * The covariance matrix adaption rate.
             */
            double ccov;

            /**
             * The step size damping.
             */
            double step_size_damping;
        };


        /**
         * Destructor.
         */
        ~Optimizer();

    private:
        /**
         * Creates a new instance of this class with the configuration supplied as argument.
         *
         * @param config The configuration.
         */
        Optimizer(const Builder &config);

        /**
         * The CMA-ES configuration.
         */
        const Builder configuration;
    };

}

#endif // ESPECIA_OPTMIZER_H
