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

    template<class Deviate, class Decompose, class Compare, class Tracer>
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
        Optimizer(Optimizer_Builder<Deviate, Decompose, Compare, Tracer> *builder) {
        }

        friend class Optimizer_Builder<Deviate, Decompose, Compare, Tracer>;
    };

    /**
     * @todo - this is work in progress
     */
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

        Optimizer_Builder &set_parent_number(unsigned parent_number = 4) {
            this->parent_number = parent_number;
            return *this;
        }

        Optimizer_Builder &set_population_size(unsigned population_size = 8) {
            this->population_size = population_size;
            return *this;
        }

        Optimizer_Builder &set_update_modulus(unsigned update_modulus = 1) {
            this->update_modulus = update_modulus;
            return *this;
        }

        Optimizer_Builder &set_accuracy_goal(double accuracy_goal = 1.0E-06) {
            this->accuracy_goal = accuracy_goal;
            return *this;
        }

        Optimizer_Builder &set_stop_generation(unsigned long stop_generation = 1000) {
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

        size_t get_n() const {
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

}

#endif //ESPECIA_OPTMIZER_H
