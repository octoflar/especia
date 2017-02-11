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

#include "optimize.h"

namespace especia {

    class Optimizer;

    /**
     * Builds a new CMA-ES optimizer.
     */
    class Optimizer_Builder {
    public:
        /**
         * The default constructor.
         */
        Optimizer_Builder();

        /**
         * The destructor.
         */
        ~Optimizer_Builder();

        /**
         * Configures the problem dimension.
         *
         * @param n The problem dimension.
         * @return this builder.
         */
        Optimizer_Builder &with_problem_dimension(unsigned n = 1);

        /**
         * Configures the parent number.
         * @param parent_number The parent number.
         * @return this builder.
         */
        Optimizer_Builder &with_parent_number(unsigned parent_number = 4);

        /**
         * Configures the population size.
         *
         * @param population_size The population size.
         * @return this builder.
         */
        Optimizer_Builder &with_population_size(unsigned population_size = 8);

        /**
         * Configures the covariance matrix update modulus.
         *
         * @param update_modulus The update modulus.
         * @return this builder.
         */
        Optimizer_Builder &with_update_modulus(unsigned update_modulus = 1);

        /**
         * Configures the accuracy goal.
         *
         * @param accuracy_goal The accuracy goal.
         * @return this builder.
         */
        Optimizer_Builder &with_accuracy_goal(double accuracy_goal = 1.0E-04);

        /**
         * Configures the stop generation.
         *
         * @param stop_generation The stop generation.
         * @return this builder.
         */
        Optimizer_Builder &with_stop_generation(unsigned long stop_generation = 1000);

        /**
         * Builds a new CMA-ES optimizer.
         * @return the CMA-ES optimizer.
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
     * @todo - this is work in progress
     */
    class Optimizer {
    public:
        ~Optimizer() {

        }

    private:
        Optimizer(Optimizer_Builder builder) {
        }

        friend class Optimizer_Builder;
    };

}

#endif //ESPECIA_OPTMIZER_H
