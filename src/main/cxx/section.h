// Class for modeling spectroscopic data sections
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
         * @param n[in] The number of data points.
         */
        Section(size_t n);

        /**
         * Constructs a new instance of this class for a certain number of data points,
         * with given wavelength, flux, and uncertainty data.
         *
         * @param n[in] The number of data points.
         * @param wav[in] The wavelength data.
         * @param flx[in] The spectral flux data.
         * @param unc[in] The spectral flux uncertainty data.
         */
        Section(size_t n, const double wav[], const double flx[], const double unc[]);

        /**
         * Destructor.
         */
        ~Section();

        /**
         * Reads a data section from an input stream.
         *
         * @param is[in,out] The input stream.
         * @param a[in] The minimum wavelength to read.
         * @param b[in] The maximum wavelength to read.
         * @return the input stream.
         */
        std::istream &get(std::istream &is, double a = 0.0, double b = std::numeric_limits<double>::max());

        /**
         * Writes a data section to an output stream.
         *
         * @param is[in,out] The output stream.
         * @param a[in] The minimum wavelength to write.
         * @param b[in] The maximum wavelength to write.
         * @return the output stream.
         */
        std::ostream &put(std::ostream &os, double a = 0.0, double b = std::numeric_limits<double>::max()) const;

        /**
         * Returns the lower wavelength bound of this data section.
         *
         * @return the lower wavelength bound of this data section.
         */
        double lower_bound() const {
            return wav[0];
        }

        /**
         * Returns the central wavelength of this data section.
         *
         * @return the central wavelength of this data section.
         */
        double center() const {
            return 0.5 * (lower_bound() + upper_bound());
        }

        /**
         * Returns the upper wavelength bound of this data section.
         *
         * @return the upper wavelength bound of this data section.
         */
        double upper_bound() const {
            return (n > 1) ? wav[n - 1] : wav[0];
        }

        /**
         * Returns the width of this data section.
         *
         * @return the width of this data section.
         */
        double width() const {
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
         * Returns the value of the cost function resulting from having applied an optical
         * depth model to this section.
         *
         * @return the value of the cost function.
         */
        double cost() const;


        /**
         * Returns the value of the cost function as a function of a given optical depth
         * model.
         *
         * @tparam M The type of optical depth model.
         * @param model[in] The optical depth model.
         * @param r[in] The spectral resolution of the instrument.
         * @param m[in] The number of Legendre basis polynomials to model the background continuum.
         *
         * @return the value of the cost function.
         */
        template<class M>
        double cost(const M &model, double r, size_t m) const {
            using std::abs;
            using std::valarray;

            valarray<double> opt(n);
            valarray<double> atm(n);
            valarray<double> cat(n);
            valarray<double> cfl(n);
            valarray<double> tfl(n);
            valarray<double> fit(n);
            valarray<double> res(n);

            convolute(model, r, &opt[0], &atm[0], &cat[0]);
            continuum(m, &cat[0], &cfl[0]);

            double cost = 0.0;

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
         * @param a[in] The lower bound of the interval.
         * @param b[in] The upper bound of the interval.
         */
        void mask(double a, double b);

        /**
         * Applies an optical depth model to this section.
         *
         * @tparam M The type of optical depth model.
         * @param model[in] The optical depth model.
         * @param r[in] The spectral resolution of the instrument.
         * @param m[in] The number of Legendre basis polynomials to model the background continuum.
         *
         * @return a reference to this section.
         */
        template<class M>
        Section &apply(const M &model, double r, size_t m) {
            convolute(model, r, &opt[0], &atm[0], &cat[0]);
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
         * @param m[in] The number of Legendre basis polynomials to model the background continuum.
         * @param cat[in] The evaluated convoluted absorption term.
         * @param cfl[out] The evaluated background continuum flux.
         */
        void continuum(size_t m, const double cat[], double cfl[]) const throw(std::runtime_error);

        /**
         * Convolutes a given optical depth model with the instrumental line spread function.
         *
         * @tparam M The type of optical depth model.
         * @param model[in] The optical depth model.
         * @param r[in] The spectral resolution of the instrument.
         * @param opt[out] The evaluated optical depth.
         * @param atm[out] The evaluated absorption term.
         * @param cat[out] The evaluated convoluted absorption term.
         */
        template<class M>
        void convolute(const M &model, double r, double opt[], double atm[], double cat[]) const {
            using std::exp;
            using std::valarray;

            if (n > 2) {
                // The half width at half maximum (HWHM) of the instrumental profile.
                const double h = 0.5 * center() / (r * kilo);
                // The sample spacing.
                const double w = width() / (n - 1);
                // The Gaussian line spread function is cut at 4 HWHM where it is less than 10E-5.
                const size_t m = static_cast<size_t>(4.0 * (h / w)) + 1;

                valarray<double> p(m);
                valarray<double> q(m);

                for (size_t i = 0; i < m; ++i) {
                    primitive(i * w, h, p[i], q[i]);
                }

                for (size_t i = 0; i < n; ++i) {
                    opt[i] = model(wav[i]);
                    atm[i] = exp(-opt[i]);
                }

                // Convolution of the modelled flux with the instrumental line spread function.
                for (size_t i = 0; i < n; ++i) {
                    double a = 0.0;
                    double b = 0.0;

                    for (size_t j = 0; j + 1 < m; ++j) {
                        const size_t k = (i < j + 1) ? 0 : i - j - 1;
                        const size_t l = (i + j + 2 > n) ? n - 2 : i + j;
                        const double d = (atm[l + 1] - atm[l]) - (atm[k + 1] - atm[k]);

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
         * @param x[in] The abscissa value where to evaluate the primitive functions.
         * @param h[in] The half width at half maximum (HWHM) of the Gaussian line spread function.
         * @param p[out] The primitive function of g(x) evaluated at @ x.
         * @param q[out] The primitive function of x g(x) evaluated at @ x.
         */
        void primitive(double x, double h, double &p, double &q) const;

        /**
         * The observed wavelength data (arbitrary units).
         */
        std::valarray<double> wav;

        /**
         * The observed spectral flux data (arbitrary units).
         */
        std::valarray<double> flx;

        /**
         * The observed spectral flux uncertainty data (arbitrary units).
         */
        std::valarray<double> unc;

        /**
         * The selection mask.
         */
        std::valarray<bool> msk;

        /**
         * The evaluated optical depth model.
         */
        std::valarray<double> opt;

        /**
         * The evaluated absorption term.
         */
        std::valarray<double> atm;

        /**
         * The evaluated convoluted absorption term.
         */
        std::valarray<double> cat;

        /**
         * The evaluated background continuum flux.
         */
        std::valarray<double> cfl;

        /**
         * The evaluated spectral flux.
         */
        std::valarray<double> tfl;


        /**
         * The evaluated convoluted spectral flux.
         */
        std::valarray<double> fit;

        /**
         * The evaluated residual flux.
         */
        std::valarray<double> res;

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
     * @return a reference to the input stream.
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
     * @return a reference to the output stream.
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
     * @return a reference to the input stream.
     */
    std::istream &operator>>(std::istream &is, std::vector<Section> &sections);

    /**
    * The operator to write data sections to an output stream.
    *
    * @param os[in,out] The output stream.
    * @param section[in] The data sections.
    * @return a reference to the output stream.
    */
    std::ostream &operator<<(std::ostream &os, const std::vector<Section> &sections);
}

#endif // ESPECIA_SECTION_H
