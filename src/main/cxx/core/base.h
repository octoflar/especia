/// @file base.h
/// Base types, mathematical and physical constants and functions.
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
#ifndef ESPECIA_BASE_H
#define ESPECIA_BASE_H

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace especia {

    /**
     * The class of continous univariate functions @c f(x) whose derivative exists and is continous.
     */
    template<class T>
    class C1 {
    public:
        /**
         * The first argument is the abscissa @c x, the second argument returns the value of @c f(x),
         * and the third argument returns the value of the first derivative @c f'(x).
         */
        typedef void (&type)(const T &, T &, T &);

    private:
        /**
         * Private constructor prevents instantiation.
         */
        C1() {

        }
    };

    /**
     * The type of natural numbers (8 decimal digits required) including zero (denoted in maths as set N).
     */
    typedef unsigned long int L_type;

    /**
     * The type of natural numbers (4 decimal digits required) including zero (denoted in maths as set N).
     */
    typedef unsigned int N_type;

    /**
     * The type of real numbers (denoted in maths as set R).
     */
    typedef double R_type;

    /**
     * The type of natural numbers (32 binary digits required) including zero (denoted in maths as set N).
     */
    typedef unsigned long int W_type;

    /**
     * The type of integer numbers (4 decimal digits required) including zero (denoted in maths as set Z).
     */
    typedef int Z_type;

    /**
     * Pi. <https://www.wolframalpha.com/input/?i=pi+to+42+digits>
     */
    const R_type pi = 3.14159265358979323846264338327950288419717;

    /**
     * The square root of Pi. <https://www.wolframalpha.com/input/?i=square+root+of+pi+to+42+digits>
     */
    const R_type sqrt_of_pi = 1.77245385090551602729816748334114518279755;

    /**
     * The electric constant (F m-1). *NIST SP 961 (Sept/2015)*
     */
    const R_type electric_constant = 8.854187817E-12;

    /**
     * The electron mass (kg). *NIST SP 961 (Sept/2015)*
     */
    const R_type electron_mass = 9.10938356E-31;

    /**
     * The elementary charge (C). *NIST SP 961 (Sept/2015)*
     */
    const R_type elementary_charge = 1.6021766208E-19;

    /**
     * SI prefix. The spectral resolution of an instrument is expressed in units of this number.
     */
    const R_type kilo = 1.0E+03;

    /**
     * SI prefix. Variation of the fine-structure constant is expressed in units of this number.
     */
    const R_type micro = 1.0E-06;

    /**
     * The speed of light in vacuum (m s-1). *NIST SP 961 (Sept/2015)*
     */
    const R_type speed_of_light = 299792458.0;

    /**
     * Converts a numerical string into a number.
     *
     * @tparam T The number type.
     *
     * @param s The string.
     * @return the number.
     *
     * @throw invalid_argument when the string cannot be converted into a number of requested type.
     */
    template<class T>
    static T convert(const std::string &s) throw(std::invalid_argument) {
        using std::invalid_argument;
        using std::istringstream;

        istringstream iss(s);
        T t;

        if (iss >> t) {
            return t;
        }

        throw invalid_argument(
                "especia::convert(): Error: the expression '" + s + "' cannot be converted into a number");
    }

    /**
     * Returns the photon redshift as a function of relative radial velocity between observer and emitter.
     *
     * @param[in] v The relative radial velocity between observer and emitter (m s-1).
     * @return the photon redshift.
     */
    template<class T>
    static T redshift(const T &v) {
        return std::sqrt((T(1) + v / T(299792458)) / (T(1) - v / T(299792458))) - T(1);
    }

    /**
     * Returns the square of a number.
     *
     * @tparam T The number type.
     *
     * @param[in] x The number.
     * @return the square of the number.
     */
    template<class T>
    T sqr(const T &x) {
        return (x == T(0)) ? T(0) : x * x;
    }

}

#endif // ESPECIA_BASE_H
