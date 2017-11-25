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
        explicit Section(size_t n_in);

        /**
         * Constructs a new instance of this class for a certain number of data points,
         * with given wavelength, flux, and uncertainty data.
         *
         * @param[in] n_in The number of data points.
         * @param[in] wav The wavelength data.
         * @param[in] flx The spectral flux data.
         * @param[in] unc The spectral flux uncertainty data.
         */
        Section(size_t n_in, const real wav[], const real flx[], const real unc[]);

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
        std::istream &get(std::istream &is, real a = 0.0, real b = std::numeric_limits<real>::max());

        /**
         * Writes a data section to an output stream.
         *
         * @param[in,out] is The output stream.
         * @param[in] a The minimum wavelength to write.
         * @param[in] b The maximum wavelength to write.
         * @return the output stream.
         */
        std::ostream &put(std::ostream &os, real a = 0.0, real b = std::numeric_limits<real>::max()) const;

        /**
         * Returns the lower wavelength bound of this data section.
         *
         * @return the lower wavelength bound of this data section.
         */
        real lower_bound() const {
            return n > 0 ? wav[0] : 0.0;
        }

        /**
         * Returns the central wavelength of this data section.
         *
         * @return the central wavelength of this data section.
         */
        real center() const {
            return 0.5 * (lower_bound() + upper_bound());
        }

        /**
         * Returns the upper wavelength bound of this data section.
         *
         * @return the upper wavelength bound of this data section.
         */
        real upper_bound() const {
            return (n > 0) ? wav[n - 1] : 0.0;
        }

        /**
         * Returns the width of this data section.
         *
         * @return the width of this data section.
         */
        real width() const {
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
        real cost() const;


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
        real cost(const T &tau, real r, natural m) const {
            using std::abs;
            using std::valarray;

            valarray<real> opt(n);
            valarray<real> atm(n);
            valarray<real> cat(n);
            valarray<real> cfl(n);
            valarray<real> tfl(n);
            valarray<real> fit(n);
            valarray<real> res(n);

            convolute(r, tau, opt, atm, cat);
            continuum(m, cat, cfl);

            tfl = cfl * atm;
            fit = cfl * cat;
            res = (flx - fit) / unc;

            real cost = 0.0;
            // @todo - vectorize
            for (size_t i = 0; i < n; ++i) {
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
        void mask(real a, real b);

        /**
         * Applies an optical depth and background continuum model to this section.
         *
         * @tparam T The type of optical depth model.
         *
         * @param[in] m The number of Legendre basis polynomials to model the background continuum.
         * @param[in] r The spectral resolution of the instrument.
         * @param[in] tau The optical depth model.
         *
         * @return this section.
         */
        template<class T>
        Section &apply(natural m, real r, const T &tau) {
            convolute(r, tau, opt, atm, cat);
            continuum(m, cat, cfl);

            tfl = cfl * atm;
            fit = cfl * cat;
            res = (flx - fit) / unc;

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
        void continuum(natural m, const std::valarray<real> &cat, std::valarray<real> &cfl) const;

        /**
         * Convolutes a given optical depth model with the instrumental line spread function.
         *
         * @tparam T The type of optical depth model.
         *
         * @param[in] r The spectral resolution of the instrument.
         * @param[in] tau The optical depth model.
         * @param[out] opt The evaluated optical depth.
         * @param[out] atm The evaluated absorption term.
         * @param[out] cat The evaluated convoluted absorption term.
         */
        template<class T>
        void convolute(real r, const T &tau, std::valarray<real> &opt, std::valarray<real> &atm,
                       std::valarray<real> &cat) const {
            using std::exp;
            using std::valarray;

            if (n > 2) {
                // The half width at half maximum (HWHM) of the instrumental profile.
                const real h = 0.5 * center() / (r * kilo);
                // The sample spacing.
                const real w = width() / (n - 1);
                // The Gaussian line spread function is truncated at 4 HWHM where it is less than 10E-5.
                const natural m = static_cast<natural>(4.0 * (h / w)) + 1;

                valarray<real> p(m);
                valarray<real> q(m);

                for (natural i = 0; i < m; ++i) {
                    // @todo - vectorize
                    primitive(i * w, h, p[i], q[i]);
                }

                for (size_t i = 0; i < n; ++i) {
                    // @todo - vectorize
                    opt[i] = tau(wav[i]);
                }
                atm = exp(-opt);

                // Convolution of the modelled flux with the instrumental line spread function.
                for (size_t i = 0; i < n; ++i) {
                    real a = 0.0;
                    real b = 0.0;

                    // @todo - vectorize, if possible
                    for (natural j = 0; j + 1 < m; ++j) {
                        const size_t k = (i < j + 1) ? 0 : i - j - 1;
                        const size_t l = (i + j + 2 > n) ? n - 2 : i + j;
                        const real d = (atm[l + 1] - atm[l]) - (atm[k + 1] - atm[k]);

                        a += (p[j + 1] - p[j]) * (atm[k + 1] + atm[l] - real(j) * d);
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
        void primitive(real x, real h, real &p, real &q) const;

        /**
         * The observed wavelength data (arbitrary units).
         */
        std::valarray<real> wav;

        /**
         * The observed spectral flux data (arbitrary units).
         */
        std::valarray<real> flx;

        /**
         * The observed spectral flux uncertainty data (arbitrary units).
         */
        std::valarray<real> unc;

        /**
         * The selection mask.
         */
        std::valarray<bool> msk;

        /**
         * The evaluated optical depth model.
         */
        std::valarray<real> opt;

        /**
         * The evaluated absorption term.
         */
        std::valarray<real> atm;

        /**
         * The evaluated convoluted absorption term.
         */
        std::valarray<real> cat;

        /**
         * The evaluated background continuum flux.
         */
        std::valarray<real> cfl;

        /**
         * The evaluated spectral flux.
         */
        std::valarray<real> tfl;

        /**
         * The evaluated convoluted spectral flux.
         */
        std::valarray<real> fit;

        /**
         * The evaluated residual flux.
         */
        std::valarray<real> res;

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

