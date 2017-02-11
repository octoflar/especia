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

    class Optimizer_Builder;

    /**
     * @todo - this is work in progress
     */
    template<class Deviate, class Decompose, class Compare, class Tracer>
    class Optimizer {
    public:
        ~Optimizer() {

        }

    private:
        Optimizer(Deviate &deviate, Decompose &decompose, Compare &compare, Tracer &tracer,
                  Optimizer_Builder &builder) {
        }

        friend class Optimizer_Builder;
    };

    /**
     * @todo - this is work in progress
     */
    class Optimizer_Builder {
    public:
        Optimizer_Builder();

        ~Optimizer_Builder();

        Optimizer_Builder &with_problem_dimension(unsigned n = 1) {
            this->n = n;
            set_strategy_parameters();
            return *this;
        }

        Optimizer_Builder &with_parent_number(unsigned parent_number = 4);

        Optimizer_Builder &with_population_size(unsigned population_size = 8);

        Optimizer_Builder &with_update_modulus(unsigned update_modulus = 1) {
            this->update_modulus = update_modulus;
            return *this;
        }

        Optimizer_Builder &with_accuracy_goal(double accuracy_goal = 1.0E-06) {
            this->accuracy_goal = accuracy_goal;
            return *this;
        }

        Optimizer_Builder &with_stop_generation(unsigned long stop_generation = 1000) {
            this->stop_generation = stop_generation;
            return *this;
        }

        template<class Deviate, class Decompose, class Compare, class Tracer>
        Optimizer<Deviate, Decompose, Compare, Tracer>
        build(Deviate &deviate, Decompose &decompose, Compare &compare, Tracer &tracer) {
            return Optimizer<Deviate, Decompose, Compare, Tracer>(deviate, decompose, compare, tracer, *this);
        };

        size_t get_problem_dimension() const {
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

        std::valarray<double> get_weights() const {
            return w;
        }

        double get_step_size_cumulation_rate() const {
            return cs;
        }

        double get_distribution_cumulation_rate() const {
            return cc;
        }

        double get_covariance_matrix_adaption_rate() const {
            return ccov;
        }

        double get_covariance_matrix_adaption_mixing() const {
            return acov;
        }

        double get_step_size_damping() const {
            return step_size_damping;
        }

    private:

        void set_strategy_parameters();

        size_t n;

        unsigned int parent_number;
        unsigned int population_size;
        unsigned int update_modulus;

        double accuracy_goal;
        unsigned long stop_generation;

        std::valarray<double> w;

        double wv;
        double cs;
        double cc;
        double acov;
        double ccov;
        double step_size_damping;

    };

}

#endif //ESPECIA_OPTMIZER_H
