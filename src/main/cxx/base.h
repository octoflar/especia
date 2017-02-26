// Basic constants and functions
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

#include <valarray>

namespace especia {

    /**
     * Pi. <https://www.wolframalpha.com/input/?i=pi+to+42+digits>
     */
    const double pi = 3.14159265358979323846264338327950288419717;

    /**
     * The square root of Pi. <https://www.wolframalpha.com/input/?i=square+root+of+pi+to+42+digits>
     */
    const double sqrt_of_pi = 1.77245385090551602729816748334114518279755;

    /**
     * The electric constant (F m-1). *NIST SP 961 (Sept/2015)*
     */
    const double electric_constant = 8.854187817E-12;

    /**
     * The electron mass (kg). *NIST SP 961 (Sept/2015)*
     */
    const double electron_mass = 9.10938356E-31;

    /**
     * The elementary charge (C). *NIST SP 961 (Sept/2015)*
     */
    const double elementary_charge = 1.6021766208E-19;

    /**
     * SI prefix. The spectral resolution of an instrument is expressed in units of this number.
     */
    const double kilo = 1.0E+03;

    /**
     * SI prefix. Variation of the fine-structure constant is expressed in units of this number.
     */
    const double micro = 1.0E-06;

    /**
     * The speed of light in vacuum (m s-1). *NIST SP 961 (Sept/2015)*
     */
    const double speed_of_light = 299792458.0;

    /**
     * Returns the square of a number.
     *
     * @tparam Number The number type.
     *
     * @param x[in] The number.
     * @return the square of the number.
     */
    template<class Number>
    inline
    Number sqr(const Number &x) {
        return (x == Number(0)) ? Number(0) : x * x;
    }

    /**
     * Returns the relativistic Doppler factor for a given radial velocity.
     *
     * @param v[in] The radial velocity (m s-1).
     * @return the Doppler factor.
     */
    inline
    double dop(const double &v) {
        return sqrt((1.0 + v / speed_of_light) / (1.0 - v / speed_of_light));
    }

    /**
     * Used to convert wavelength in vacuum to wavelength in air.
     *
     * Further reading:
     *
     * B. Edlen (1966).
     *   *The refractive index of air.*
     *   Metrologia, 2, 2, 71-80.
     *   http://dx.doi.org/10.1088/0026-1394/2/2/002
     *
     * @param x[in] The quotient 10 / wavelength in vacuum.
     * @return the quotient 10 / wavelength in air.
     */
    inline
    double edlen(const double &x) {
        return (1.0000834213 + 1.5997E-10 / (0.0000389 - x * x) + 2.406030E-08 / (0.000130 - x * x)) * x;
    }

    /**
     * Used to convert wavelength in air to wavelength in vacuum.
     *
     * Further reading:
     *
     * B. Edlen (1966).
     *   *The refractive index of air.*
     *   Metrologia, 2, 2, 71-80.
     *   http://dx.doi.org/10.1088/0026-1394/2/2/002
     *
     * @param x[in] The quotient 10 / wavelength in vacuum.
     * @param y[out] The quotient 10 / wavelength in air.
     * @param z[out] The derivative of @c y with respect to @c x.
     */
    inline
    void edlen(const double &x, double &y, double &z) {
        const double a = 1.0000834213 + 1.5997E-10 / (0.0000389 - x * x) + 2.406030E-08 / (0.000130 - x * x);
        const double b = x * ((3.1994E-10 * x) / sqr(0.0000389 - x * x) + (4.81206E-08 * x) / sqr(0.000130 - x * x));

        y = a * x;
        z = a + b;
    }

    /**
     * Used to convert wavelength in vacuum to wavelength in air.
     *
     * Further reading:
     *
     * B. Edlen (1953).
     *   *The dispersion of standard air.*
     *   Journal of the Optical Society of America, 43, 5, 339.
     *
     * @param x[in] The quotient 10 / wavelength in vacuum.
     * @return the quotient 10 / wavelength in air.
     */
    inline
    double edlen_1953(const double &x) {
        return (1.0000643280 + 2.5540E-10 / (0.0000410 - x * x) + 2.949810E-08 / (0.000146 - x * x)) * x;
    }

    /**
     * Used to convert wavelength in air to wavelength in vacuum.
     *
     * Further reading:
     *
     * B. Edlen (1953).
     *   *The dispersion of standard air.*
     *   Journal of the Optical Society of America, 43, 5, 339.
     *
     * @param x[in] The quotient 10 / wavelength in vacuum.
     * @param y[out] The quotient 10 / wavelength in air.
     * @param z[out] The derivative of @c y with respect to @c x.
     */
    inline
    void edlen_1953(const double &x, double &y, double &z) {
        const double a = 1.0000643280 + 2.5540E-10 / (0.0000410 - x * x) + 2.949810E-08 / (0.000146 - x * x);
        const double b = x * ((5.1080E-10 * x) / sqr(0.0000410 - x * x) + (5.89962E-08 * x) / sqr(0.000146 - x * x));

        y = a * x;
        z = a + b;
    }

    /**
     * An indirect comparing.
     *
     * @tparam Number The base value type.
     * @tparam Compare The strategy to compare base values directly.
     */
    template<class Number, class Compare>
    class Indirect_Compare {
    public:
        /**
         * Constructs a new indirect comparing.
         *
         * @param[in] v The base values.
         * @param[in] c The direct base value comparing.
         */
        Indirect_Compare(const std::valarray<Number> &v, const Compare &c)
                : values(v), compare(c) {
        }

        /**
         * Destructor.
         */
        ~Indirect_Compare() {
        }

        /**
         * The indirect comparing operator.
         *
         * @param[in] i An index into the set of base values.
         * @param[in] j An index into the set of base values.
         * @return the result of comparing the indexed base values directly.
         */
        bool operator()(const size_t &i, const size_t &j) const {
            return compare(values[i], values[j]);
        }

    private:
        const std::valarray<Number> &values;
        const Compare &compare;
    };

}

#endif // ESPECIA_BASE_H
