/// @file section.h
/// Class for modeling spectroscopic data sections.
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
#ifndef ESPECIA_SECTION_H
#define ESPECIA_SECTION_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <valarray>
#include <vector>

#include "base.h"

namespace especia {

    /**
     * Represents a section of (observed and modelled) spectroscopic data.
     */
    class Section {
    public:

        /**
         * Constructs a new instance of this class, which contains no data.
         */
        Section();

        /**
         * Constructs a new instance of this class for a certain number of data points.
         *
         * @param[in] n_in The number of data points.
         */
        Section(size_t n_in);

        /**
         * Constructs a new instance of this class for a certain number of data points,
         * with given wavelength, flux, and uncertainty data.
         *
         * @param[in] n_in The number of data points.
         * @param[in] wav The wavelength data.
         * @param[in] flx The spectral flux data.
         * @param[in] unc The spectral flux uncertainty data.
         */
        Section(size_t n_in, const R_type wav[], const R_type flx[], const R_type unc[]);

        /**
         * The destructor.
         */
        ~Section();

        /**
         * Reads a data section from an input stream.
         *
         * @param[in,out] is The input stream.
         * @param[in] a The minimum wavelength to read.
         * @param[in] b The maximum wavelength to read.
         * @return the input stream.
         */
        std::istream &get(std::istream &is, R_type a = 0.0, R_type b = std::numeric_limits<R_type>::max());

        /**
         * Writes a data section to an output stream.
         *
         * @param[in,out] is The output stream.
         * @param[in] a The minimum wavelength to write.
         * @param[in] b The maximum wavelength to write.
         * @return the output stream.
         */
        std::ostream &put(std::ostream &os, R_type a = 0.0, R_type b = std::numeric_limits<R_type>::max()) const;

        /**
         * Returns the lower wavelength bound of this data section.
         *
         * @return the lower wavelength bound of this data section.
         */
        R_type lower_bound() const {
            return n > 0 ? wav[0] : 0.0;
        }

        /**
         * Returns the central wavelength of this data section.
         *
         * @return the central wavelength of this data section.
         */
        R_type center() const {
            return 0.5 * (lower_bound() + upper_bound());
        }

        /**
         * Returns the upper wavelength bound of this data section.
         *
         * @return the upper wavelength bound of this data section.
         */
        R_type upper_bound() const {
            return (n > 0) ? wav[n - 1] : 0.0;
        }

        /**
         * Returns the width of this data section.
         *
         * @return the width of this data section.
         */
        R_type width() const {
            return upper_bound() - lower_bound();
        }

        /**
         * Returns the number of data points in this section.
         *
         * @return the number of data points.
         */
        size_t data_count() const {
            return n;
        }

        /**
         * Returns the number of valid data points in this section.
         *
         * @return the number of valid data points.
         */
        size_t valid_data_count() const;

        /**
         * Returns the current value of the cost function.
         *
         * @return the current value of the cost function.
         */
        R_type cost() const;


        /**
         * Returns the value of the cost function as a function of a given optical depth
         * model.
         *
         * @tparam T The type of optical depth model.
         *
         * @param[in] tau The optical depth model.
         * @param[in] r The spectral resolution of the instrument.
         * @param[in] m The number of Legendre basis polynomials to model the background continuum.
         *
         * @return the value of the cost function.
         *
         * @remark calling this method is thread safe, if the optical depth model is thread safe.
         */
        template<class T>
        R_type cost(const T &tau, R_type r, N_type m) const {
            using std::abs;
            using std::valarray;

            valarray<R_type> opt(n);
            valarray<R_type> atm(n);
            valarray<R_type> cat(n);
            valarray<R_type> cfl(n);
            valarray<R_type> tfl(n);
            valarray<R_type> fit(n);
            valarray<R_type> res(n);

            convolute(tau, r, &opt[0], &atm[0], &cat[0]);
            continuum(m, &cat[0], &cfl[0]);

            R_type cost = 0.0;

            for (size_t i = 0; i < n; ++i) {
                tfl[i] = cfl[i] * atm[i];
                fit[i] = cfl[i] * cat[i];

                res[i] = (flx[i] - fit[i]) / unc[i];

                if (msk[i]) {
                    cost += res[i] * res[i];
                }
            }

            return 0.5 * cost;
        }

        /**
         * Masks the data in a certain interval as invalid.
         *
         * @param[in] a The lower bound of the interval.
         * @param[in] b The upper bound of the interval.
         */
        void mask(R_type a, R_type b);

        /**
         * Applies an optical depth model to this section.
         *
         * @tparam T The type of optical depth model.
         *
         * @param[in] tau The optical depth model.
         * @param[in] r The spectral resolution of the instrument.
         * @param[in] m The number of Legendre basis polynomials to model the background continuum.
         *
         * @return this section.
         */
        template<class T>
        Section &apply(const T &tau, R_type r, N_type m) {
            convolute(tau, r, &opt[0], &atm[0], &cat[0]);
            continuum(m, &cat[0], &cfl[0]);

            for (size_t i = 0; i < n; ++i) {
                tfl[i] = cfl[i] * atm[i];
                fit[i] = cfl[i] * cat[i];

                res[i] = (flx[i] - fit[i]) / unc[i];
            }
            return *this;
        }

    private:
        /**
         * Calculates an optimized background continuum.
         *
         * @param[in] m The number of Legendre basis polynomials to model the background continuum.
         * @param[in] cat The evaluated convoluted absorption term.
         * @param[out] cfl The evaluated background continuum flux.
         */
        void continuum(N_type m, const R_type cat[], R_type cfl[]) const throw(std::runtime_error);

        /**
         * Convolutes a given optical depth model with the instrumental line spread function.
         *
         * @tparam T The type of optical depth model.
         *
         * @param[in] tau The optical depth model.
         * @param[in] r The spectral resolution of the instrument.
         * @param[out] opt The evaluated optical depth.
         * @param[out] atm The evaluated absorption term.
         * @param[out] cat The evaluated convoluted absorption term.
         */
        template<class T>
        void convolute(const T &tau, R_type r, R_type opt[], R_type atm[], R_type cat[]) const {
            using std::exp;
            using std::valarray;

            if (n > 2) {
                // The half width at half maximum (HWHM) of the instrumental profile.
                const R_type h = 0.5 * center() / (r * kilo);
                // The sample spacing.
                const R_type w = width() / (n - 1);
                // The Gaussian line spread function is truncated at 4 HWHM where it is less than 10E-5.
                const N_type m = static_cast<N_type>(4.0 * (h / w)) + 1;

                valarray<R_type> p(m);
                valarray<R_type> q(m);

                for (N_type i = 0; i < m; ++i) {
                    primitive(i * w, h, p[i], q[i]);
                }

                for (size_t i = 0; i < n; ++i) {
                    opt[i] = tau(wav[i]);
                    atm[i] = exp(-opt[i]);
                }

                // Convolution of the modelled flux with the instrumental line spread function.
                for (size_t i = 0; i < n; ++i) {
                    R_type a = 0.0;
                    R_type b = 0.0;

                    for (N_type j = 0; j + 1 < m; ++j) {
                        const size_t k = (i < j + 1) ? 0 : i - j - 1;
                        const size_t l = (i + j + 2 > n) ? n - 2 : i + j;
                        const R_type d = (atm[l + 1] - atm[l]) - (atm[k + 1] - atm[k]);

                        a += (p[j + 1] - p[j]) * (atm[k + 1] + atm[l] - j * d);
                        b += (q[j + 1] - q[j]) * d;
                    }

                    cat[i] = a + b / w;
                }
            }
        }

        /**
         * Evaluates the primitive functions of g(x) and x g(x), where g(x) is the (Gaussian)
         * line spread function of the instrument.
         *
         * @param[in] x The abscissa value where to evaluate the primitive functions.
         * @param[in] h The half width at half maximum (HWHM) of the Gaussian line spread function.
         * @param[out] p The primitive function of g(x) evaluated at @ x.
         * @param[out] q The primitive function of x g(x) evaluated at @ x.
         */
        void primitive(R_type x, R_type h, R_type &p, R_type &q) const;

        /**
         * The observed wavelength data (arbitrary units).
         */
        std::valarray<R_type> wav;

        /**
         * The observed spectral flux data (arbitrary units).
         */
        std::valarray<R_type> flx;

        /**
         * The observed spectral flux uncertainty data (arbitrary units).
         */
        std::valarray<R_type> unc;

        /**
         * The selection mask.
         */
        std::valarray<bool> msk;

        /**
         * The evaluated optical depth model.
         */
        std::valarray<R_type> opt;

        /**
         * The evaluated absorption term.
         */
        std::valarray<R_type> atm;

        /**
         * The evaluated convoluted absorption term.
         */
        std::valarray<R_type> cat;

        /**
         * The evaluated background continuum flux.
         */
        std::valarray<R_type> cfl;

        /**
         * The evaluated spectral flux.
         */
        std::valarray<R_type> tfl;

        /**
         * The evaluated convoluted spectral flux.
         */
        std::valarray<R_type> fit;

        /**
         * The evaluated residual flux.
         */
        std::valarray<R_type> res;

        /**
         * The number of data points.
         */
        size_t n;
    };


    /**
     * The operator to read a data section from an input stream.
     *
     * @param is[in,out] The input stream.
     * @param section[out] The data section.
     * @return the input stream.
     */
    inline
    std::istream &operator>>(std::istream &is, Section &section) {
        return section.get(is);
    }

    /**
     * The operator to write a data section to an output stream.
     *
     * @param os[in,out] The output stream.
     * @param section[in] The data section.
     * @return the output stream.
     */
    inline
    std::ostream &operator<<(std::ostream &os, const Section &section) {
        return section.put(os);
    }

    /**
     * The operator to read data sections from an input stream.
     *
     * @param is[in,out] The input stream.
     * @param section[out] The data sections.
     * @return the input stream.
     */
    std::istream &operator>>(std::istream &is, std::vector<Section> &sections);

    /**
    * The operator to write data sections to an output stream.
    *
    * @param os[in,out] The output stream.
    * @param section[in] The data sections.
    * @return the output stream.
    */
    std::ostream &operator<<(std::ostream &os, const std::vector<Section> &sections);
}

#endif // ESPECIA_SECTION_H
