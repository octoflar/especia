// Class for modeling absorption line regions
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
    // Class for modeling absorption line regions
    class Section;

    std::istream &operator>>(std::istream &is, Section &s);

    std::istream &operator>>(std::istream &is, std::vector<Section> &s);

    std::ostream &operator<<(std::ostream &os, const Section &s);

    std::ostream &operator<<(std::ostream &os, const std::vector<Section> &s);
}

class especia::Section {
public:
    Section();

    Section(size_t n);

    Section(const double wav[], const double flx[], const double err[], size_t n);

    ~Section();

    std::istream &get(std::istream &is, double a = 0.0, double b = std::numeric_limits<double>::max());

    std::ostream &put(std::ostream &os, double a = 0.0, double b = std::numeric_limits<double>::max()) const;

    double begin() const;

    double center() const;

    double end() const;

    double width() const;

    size_t size() const;

    size_t selection_size() const;

    double cost() const;

    template<class optical_depth>
    double cost(const optical_depth &t, size_t m, double r) const;

    void mask(double a, double b);

    template<class optical_depth>
    void compute_model(const optical_depth &t, size_t m, double r);

private:
    template<class optical_depth>
    void convolve(const optical_depth &t, double r, double *opt, double *atm, double *cat) const;

    void integrals(double x, double fwhm, double &p, double &q) const;
        // the indefinite integrals of P(x) and xP(x), where P(x) is the instrumental profile

    void continuum(size_t m, const double cat[], double cfl[]) const throw(std::runtime_error);

    std::valarray<double> wav; // wavelength data
    std::valarray<double> flx; // flux data
    std::valarray<double> err; // flux data uncertainty

    std::valarray<bool> msk; // selection mask

    std::valarray<double> opt; // true optical depth
    std::valarray<double> atm; // absorption term
    std::valarray<double> cat; // convoluted absorption term
    std::valarray<double> cfl; // continuum flux
    std::valarray<double> tfl; // true flux
    std::valarray<double> fit; // convoluted true flux
    std::valarray<double> res; // residuals

    size_t n; // number of data points
};

inline
double especia::Section::begin() const {
    return wav[0];
}

inline
double especia::Section::center() const {
    return 0.5 * (begin() + end());
}

inline
double especia::Section::end() const {
    return (n > 1) ? wav[n - 1] : wav[0];
}

inline
double especia::Section::width() const {
    return end() - begin();
}

inline
size_t especia::Section::size() const {
    return n;
}

template<class optical_depth>
void especia::Section::convolve(const optical_depth &t, double r, double *opt, double *atm, double *cat) const {
    using std::exp;
    using std::valarray;

    if (n > 2) {
        const double hwhm = 0.5 * center() / (r * kilo);
            // HWHM (half width at half maximum) of the instrumental profile
        const double h = width() / (n - 1);
            // sample spacing
        const size_t m = static_cast<size_t>(4.0 * (hwhm / h)) + 1;
            // cut the Gaussian profile at 4 HWHM where it is 10E-5
        valarray<double> p(m);
        valarray<double> q(m);

        for (size_t i = 0; i < m; ++i)
            integrals(i * h, hwhm, p[i], q[i]);

        for (size_t i = 0; i < n; ++i) {
            opt[i] = t(wav[i]);
            atm[i] = exp(-opt[i]);
        }

        // Convolve the true flux with the instrumental profile
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

            cat[i] = a + b / h;
        }
    }
}

template<class optical_depth>
double especia::Section::cost(const optical_depth &t, size_t m, double r) const {
    using std::abs;
    using std::valarray;

    valarray<double> opt(n);
    valarray<double> atm(n);
    valarray<double> cat(n);
    valarray<double> cfl(n);
    valarray<double> tfl(n);
    valarray<double> fit(n);
    valarray<double> res(n);

    convolve(t, r, &opt[0], &atm[0], &cat[0]);
    continuum(m, &cat[0], &cfl[0]);

    double a = 0.0;

    for (size_t i = 0; i < n; ++i) {
        tfl[i] = cfl[i] * atm[i];
        fit[i] = cfl[i] * cat[i];

        res[i] = (flx[i] - fit[i]) / err[i];

        if (msk[i])
            a += res[i] * res[i];
    }

    return 0.5 * a;
}

template<class optical_depth>
void especia::Section::compute_model(const optical_depth &t, size_t m, double r) {
    convolve(t, r, &opt[0], &atm[0], &cat[0]);
    continuum(m, &cat[0], &cfl[0]);

    for (size_t i = 0; i < n; ++i) {
        tfl[i] = cfl[i] * atm[i];
        fit[i] = cfl[i] * cat[i];

        res[i] = (flx[i] - fit[i]) / err[i];
    }
}

inline
std::istream &especia::operator>>(std::istream &is, Section &s) {
    return s.get(is);
}

inline
std::ostream &especia::operator<<(std::ostream &os, const Section &s) {
    return s.put(os);
}

#endif // ESPECIA_SECTION_H
