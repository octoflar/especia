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
     * The type of integer numbers including zero (denoted in maths as set Z).
     */
    typedef int integer;

    /**
     * The type of natural numbers including zero (denoted in maths as set N).
     */
    typedef unsigned natural;

    /**
     * The type of real numbers (denoted in maths as set R).
     */
    typedef double real;

    /**
     * The type of binary numbers (32 binary digits required).
     */
    typedef unsigned long word;

    /**
     * The class of continuous univariate functions @c f(x) whose derivative exists and is continous.
     */
    template<class T = real>
    class C1 {
    public:
        /**
         * The first argument is the abscissa @c x, the second argument returns the value of @c f(x),
         * and the third argument returns the value of the first derivative @c f'(x).
         */
        typedef void (&type)(const T &, T &, T &);

    private:
        C1() {
            // private constructor prevents instantiation
        }
    };

    /**
     * Pi. <https://www.wolframalpha.com/input/?i=pi+to+49+digits>
     */
    const real pi = real(3.141592653589793238462643383279502884197169399375L);

    /**
     * The square root of Pi. <https://www.wolframalpha.com/input/?i=square+root+of+pi+to+49+digits>
     */
    const real sqrt_of_pi = real(1.772453850905516027298167483341145182797549456123L);

    /**
     * The electric constant (F m-1). *NIST SP 961 (Sept/2015)*
     */
    const real electric_constant = 8.854187817E-12;

    /**
     * The electron mass (kg). *NIST SP 961 (Sept/2015)*
     */
    const real electron_mass = 9.10938356E-31;

    /**
     * The elementary charge (C). *NIST SP 961 (Sept/2015)*
     */
    const real elementary_charge = 1.6021766208E-19;

    /**
     * SI prefix. The spectral resolution of an instrument is expressed in units of this number.
     */
    const real kilo = 1.0E+03;

    /**
     * SI prefix. Variation of the fine-structure constant is expressed in units of this number.
     */
    const real micro = 1.0E-06;

    /**
     * The speed of light in vacuum (m s-1). *NIST SP 961 (Sept/2015)*
     */
    const real speed_of_light = 299792458.0;

    /**
     * Converts a numeric character string into a number.
     *
     * @tparam T The number type.
     *
     * @param s The string.
     * @return the number.
     *
     * @throw invalid_argument when the string cannot be converted into a number of requested type.
     */
    template<class T>
    T convert(const std::string &s) {
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
    template<class T = real>
    T redshift(const T &v) {
        return std::sqrt((T(1.0) + v / T(299792458.0)) / (T(1.0) - v / T(299792458.0))) - T(1.0);
    }

    /**
     * Solves the equation f(x) = c by means of Newton's method.
     *
     * @tparam T The number type.
     *
     * @param[in] f The function.
     * @param[in] c The constant on the right-hand side of the equation.
     * @param[in] x The initial guess of the solution.
     * @param[in] accuracy_goal The accuracy goal.
     * @param[in] max_iteration The maximum number of iterations (optional).
     *
     * @return the solution to the equation f(x) = c.
     *
     * @throw runtime_error when the accuracy goal was not reached within the prescribed number of iterations.
     */
    template<class T = real>
    T solve(typename C1<T>::type &f, T c, T x, T accuracy_goal, natural max_iteration = 100) {
        using std::abs;
        using std::runtime_error;

        T d, y, z;

        for (natural i = 0; i < max_iteration; ++i) {
            f(x, y, z);
            d = (y - c) / z;
            x -= d;
            if (abs(d) < accuracy_goal * x) {
                return x;
            }
        }

        throw runtime_error("especia::solve() Error: the required accuracy goal was not reached");
    }

    /**
     * Returns the square of a number.
     *
     * @tparam T The number type.
     *
     * @param[in] x The number.
     * @return the square of the number.
     */
    template<class T = real>
    T sq(const T &x) {
        return x * x;
    }

}

#endif // ESPECIA_BASE_H
