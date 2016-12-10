// Configuration constants
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
#ifndef ESPECIA_BASE_H
#define ESPECIA_BASE_H

namespace especia {

    /**
     * Pi.
     */
    const double pi = 3.1415926535897932384626433832795028841972; // -

    /**
     * The square root of Pi.
     */
    const double sqrt_of_pi = 1.7724538509055160272981674833411451827975; // +

    /**
     * The speed of light in vacuum (km s-1).
     */
    const double speed_of_light = 299792.458;

    /**
     * Returns the Doppler factor.
     *
     * @tparam number The number type.
     *
     * @param v The relative radial velocity between emitter and observer (km s-1)
     * @return the Doppler factor.
     */
    template<class number>
    double doppler_factor(const number &v) {
        return std::sqrt((number(1) + v / speed_of_light) / (number(1) - v / speed_of_light));
    }

    /**
     * Returns the square of a number.
     *
     * @tparam number The number type.
     *
     * @param x The number.
     * @return the square of the number.
     */
    template<class number>
    number sqr(const number &x) {
        return (x == number(0)) ? number(0) : x * x;
    }

}

#endif // ESPECIA_BASE_H
