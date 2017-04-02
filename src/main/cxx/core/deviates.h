/// @file deviates.h
/// Function-like class templates for generating various random deviates.
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
#ifndef ESPECIA_DEVIATES_H
#define ESPECIA_DEVIATES_H

#include <cmath>

#include "base.h"

namespace especia {

   /**
    * A class template to generate random normal deviates.
    *
    * The algorithm uses the polar method (e.g. Knuth, 1998,
    * Sec. 3.4.1, Algorithm P) to generate standard normally distributed random
    * deviates.
    *
    * Further reading:
    *
    * D. Knuth (1998).
    *   *The art of computer programming 2. Seminumerical algorithms.*
    *   Addison Wesley Longman, ISBN 0-201-89684-2.
    *
    * @tparam U The strategy to generate random uniform deviates.
    */
    template<class U>
    class Normal_Deviate {
    public:
        /**
         * Constructs a new instance of this class from a seed.
         *
         * @param[in] seed The seed.
         */
        Normal_Deviate(W_type seed = 5489)
                : uniform_deviate(seed), status(false) {
        }

        /**
         * Constructs a new instance of this class from a uniform deviate.
         *
         * @param[in] u The instance of this class to be copied.
         */
        Normal_Deviate(const U &u)
                : uniform_deviate(u), status(false) {
        }

        /**
         * The destructor.
         */
        ~Normal_Deviate() {
        }

        /**
         * Returns a normal random number.
         *
         * @return a normal random number.
         */
        R_type operator()() const {
            using std::log;
            using std::sqrt;

            status = !status;

            if (status) {
                R_type t;

                do {
                    x = 2.0 * uniform_deviate() - 1.0;
                    y = 2.0 * uniform_deviate() - 1.0;
                    t = x * x + y * y;
                } while (t >= 1.0 or t == 0.0);

                t = sqrt(-2.0 * (log(t) / t));
                x *= t;
                y *= t;

                return x;
            } else {
                return y;
            }
        }

        /**
         * Resets this functor with a seed.
         *
         * @param[in] seed The seed.
         */
        void reset(W_type seed = 5489) {
            uniform_deviate.reset(seed);
            status = false;
        }

        /**
         * Resets this functor with another instance.
         *
         * @param[in] u The other instance.
         */
        void reset(const U &u) {
            uniform_deviate = u;
            status = false;
        }

    private:
        const U uniform_deviate;

        mutable bool status;
        mutable R_type x, y;
    };

}

#endif // ESPECIA_DEVIATES_H
