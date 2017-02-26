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
     * Used to convert photon wavelength in vacuum to photon wavelength in air.
     *
     * Further reading:
     *
     * Donald C. Morton (2000).
     *   *Atomic Data for Resonance Absorption Lines. II. Wavelengths Longward of the Lyman Limit for Heavy Elements.*
     *   Astrophys. J. Suppl. Ser., 130, 2, 403.
     *
     * K. P. Birch and M. J. Downs (1994)
     *   *Correction to the Updated Edlén Equation for the Refractive Index of Air*
     *   Metrologia, 31, 4, 315.
     *
     * @param x[in] The wavenumber in vacuum (µm-1).
     * @return the wavenumber in air (µm-1).
     *
     * @attention The function uses wavenumber (µm-1) := 10000.0 / wavelength (Angstrom) as input and output.
     */
    inline
    double birch_1994(const double &x) {
        return (1.0 + 8.34254E-05 + 2.406147E-02 / (130.0 - x * x) + 1.5998E-04 / (38.9 - x * x)) * x;
    }

    /**
     * Used to convert photon wavelength in vacuum to photon wavelength in air (by means of Newton's method).
     *
     * Further reading:
     *
     * Donald C. Morton (2000).
     *   *Atomic Data for Resonance Absorption Lines. II. Wavelengths Longward of the Lyman Limit for Heavy Elements.*
     *   Astrophys. J. Suppl. Ser., 130, 2, 403.
     *
     * K. P. Birch and M. J. Downs (1994)
     *   *Correction to the Updated Edlén Equation for the Refractive Index of Air*
     *   Metrologia, 31, 4, 315.
     *
     * @param x[in] The wavenumber in vacuum (µm-1).
     * @param y[out] The wavenumber in air (µm-1).
     * @param z[out] The derivative of @c y with respect to @c x.
     *
     * @attention The function uses wavenumber (µm-1) := 10000.0 / wavelength (Angstrom) as input and output.
     */
    inline
    void birch_1994(const double &x, double &y, double &z) {
        const double n = 1.0 + 8.34254E-05 + 2.406147E-02 / (130.0 - x * x) + 1.5998E-04 / (38.9 - x * x);
        const double m = (4.812294E-02 * x) / sqr(130.0 - x * x) + (3.1996E-04 * x) / sqr(38.9 - x * x);

        y = x * n;
        z = n + x * m;
    }

    /**
     * Used to convert photon wavelength in vacuum to photon wavelength in air.
     *
     * Further reading:
     *
     * B. Edlén (1953).
     *   *The dispersion of standard air.*
     *   Journal of the Optical Society of America, 43, 5, 339.
     *
     * @param x[in] The wavenumber in vacuum (nm-1).
     * @return the wavenumber in air (nm-1).
     *
     * @attention The function uses wavenumber (nm-1) := 10.0 / wavelength (Angstrom) as input and output.
     *
     * @remark This formula is the IAU standard for the vacuum to standard air corrections (see resolution
     * No. C15, Commission 44, XXI General Assembly in 1991).
     */
    inline
    double edlen_1953(const double &x) {
        return (1.0000643280 + 2.5540E-10 / (0.0000410 - x * x) + 2.949810E-08 / (0.000146 - x * x)) * x;
    }

    /**
     * Used to convert photon wavelength in air to photon wavelength in vacuum (by means of Newton's method).
     *
     * Further reading:
     *
     * B. Edlén (1953).
     *   *The dispersion of standard air.*
     *   Journal of the Optical Society of America, 43, 5, 339.
     *
     * @param x[in] The wavenumber in vacuum (nm-1).
     * @param y[out] The wavenumber in air (nm-1).
     * @param z[out] The derivative of @c y with respect to @c x.
     *
     * @attention The function uses wavenumber (nm-1) := 10.0 / wavelength (Angstrom) as input and output.
     *
     * @remark This formula is the IAU standard for the vacuum to standard air corrections (see resolution
     * No. C15, Commission 44, XXI General Assembly in 1991).
     */
    inline
    void edlen_1953(const double &x, double &y, double &z) {
        const double n = 1.0000643280 + 2.5540E-10 / (0.0000410 - x * x) + 2.949810E-08 / (0.000146 - x * x);
        const double m = (5.1080E-10 * x) / sqr(0.0000410 - x * x) + (5.89962E-08 * x) / sqr(0.000146 - x * x);

        y = x * n;
        z = n + x * m;
    }

    /**
     * Used to convert photon wavelength in vacuum to photon wavelength in air.
     *
     * Further reading:
     *
     * B. Edlén (1966).
     *   *The refractive index of air.*
     *   Metrologia, 2, 2, 71-80.
     *   http://dx.doi.org/10.1088/0026-1394/2/2/002
     *
     * @param x[in] The wavenumber in vacuum (nm-1).
     * @return the wavenumber in air (nm-1).
     *
     * @attention The function uses wavenumber (nm-1) := 10.0 / wavelength (Angstrom) as input and output.
     */
    inline
    double edlen_1966(const double &x) {
        return (1.0000834213 + 1.5997E-10 / (0.0000389 - x * x) + 2.406030E-08 / (0.000130 - x * x)) * x;
    }

    /**
     * Used to convert photon wavelength in air to photon wavelength in vacuum (by means of Newton's method).
     *
     * Further reading:
     *
     * B. Edlén (1966).
     *   *The refractive index of air.*
     *   Metrologia, 2, 2, 71-80.
     *   http://dx.doi.org/10.1088/0026-1394/2/2/002
     *
     * @param x[in] The wavenumber in vacuum (nm-1).
     * @param y[out] The wavenumber in air (nm-1).
     * @param z[out] The derivative of @c y with respect to @c x.
     *
     * @attention The function uses wavenumber (nm-1) := 10.0 / wavelength (Angstrom) as input and output.
     */
    inline
    void edlen_1966(const double &x, double &y, double &z) {
        const double n = 1.0000834213 + 1.5997E-10 / (0.0000389 - x * x) + 2.406030E-08 / (0.000130 - x * x);
        const double m = (3.1994E-10 * x) / sqr(0.0000389 - x * x) + (4.81206E-08 * x) / sqr(0.000130 - x * x);

        y = x * n;
        z = n + x * m;
    }

    /**
     * An indirect comparing of indexes.
     *
     * @tparam Number The base value type.
     * @tparam Compare The strategy to compare base values directly.
     */
    template<class Number, class Compare>
    class Index_Compare {
    public:
        /**
         * Constructs a new index comparing.
         *
         * @param[in] v The base values.
         * @param[in] c The direct base value comparing.
         */
        Index_Compare(const std::valarray<Number> &v, const Compare &c)
                : values(v), compare(c) {
        }

        /**
         * Destructor.
         */
        ~Index_Compare() {
        }

        /**
         * The index comparing operator.
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
