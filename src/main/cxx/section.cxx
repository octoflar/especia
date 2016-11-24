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
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include "section.h"

RQ::section::section()
        : wav(),
          flx(),
          err(),
          msk(),
          opt(),
          atm(),
          cat(),
          cfl(),
          tfl(),
          fit(),
          res(),
          n(0) {
}

RQ::section::section(size_t m)
        : wav(0.0, m),
          flx(0.0, m),
          err(0.0, m),
          msk(1, m),
          opt(0.0, m),
          atm(0.0, m),
          cat(0.0, m),
          cfl(0.0, m),
          tfl(0.0, m),
          fit(0.0, m),
          res(0.0, m),
          n(m) {
}

RQ::section::section(const double x[], const double y[], const double z[], size_t m)
        : wav(x, m),
          flx(y, m),
          err(z, m),
          msk(1, m),
          opt(0.0, m),
          atm(0.0, m),
          cat(0.0, m),
          cfl(0.0, m),
          tfl(0.0, m),
          fit(0.0, m),
          res(0.0, m),
          n(m) {
}

RQ::section::~section() {
}

void
RQ::section::continuum(size_t m, const double cat[], double cfl[]) const throw(std::runtime_error) {
    using std::fill;
    using std::sqrt;
    using std::runtime_error;
    using std::valarray;

    if (m > 0) {
        valarray<double> b(0.0, m);
        valarray<double> c(0.0, m);
        valarray<double> l(1.0, m);

        valarray<valarray<double> > a(b, m);

        // Linear optimization problem, establish the normal equations
        for (size_t i = 0; i < n; ++i)
            if (msk[i]) {
                const double x = 2.0 * (wav[i] - wav[0]) / length() - 1.0;
                // map wavelength domain onto the interval [-1, 1]

                double l1 = 1.0;
                double l2 = 0.0;

                // Compute the higher order Legendre polynomials
                for (size_t j = 1; j < m; ++j) {
                    const double l3 = l2;

                    l2 = l1;
                    l1 = ((2 * j - 1) * x * l2 - (j - 1) * l3) / j;
                    l[j] = l1;
                }

                for (size_t j = 0; j < m; ++j) {
                    for (size_t k = j; k < m; ++k)
                        a[j][k] += (cat[i] * cat[i] * l[j] * l[k]) / (err[i] * err[i]);

                    b[j] += (flx[i] * cat[i] * l[j]) / (err[i] * err[i]);
                }
            }

        // Solve the normal equations using Cholesky decomposition (e.g. Press et al. 2002)
        for (size_t i = 0; i < m; ++i)
            for (size_t j = i; j < m; ++j) {
                double s = a[i][j];

                for (size_t k = 0; k < i; ++k)
                    s -= a[i][k] * a[j][k];

                if (i < j)
                    a[j][i] = s / a[i][i];
                else if (s > 0.0)
                    a[i][i] = sqrt(s);
                else // the normal equations are (numerically) singular
                    throw runtime_error("RQ::section::continuum(): Error: normal equations are numerically singular");
            }
        for (size_t i = 0; i < m; ++i) {
            double s = b[i];

            for (size_t k = 0; k < i; ++k)
                s -= a[i][k] * c[k];

            c[i] = s / a[i][i];
        }
        for (size_t i = m - 1; i + 1 > 0; --i) {
            double s = c[i];

            for (size_t k = i + 1; k < m; ++k)
                s -= a[k][i] * c[k];

            c[i] = s / a[i][i];
        }

        // Compute the continuum flux
        for (size_t i = 0; i < n; ++i) {
            const double x = 2.0 * (wav[i] - wav[0]) / length() - 1.0;
                // map wavelength domain onto the interval [-1, 1]

            double l1 = 1.0;
            double l2 = 0.0;

            cfl[i] = c[0];

            // Accumulate the higher order Legendre polynomials
            for (size_t k = 1; k < m; ++k) {
                const double l3 = l2;

                l2 = l1;
                l1 = ((2 * k - 1) * x * l2 - (k - 1) * l3) / k;
                cfl[i] += c[k] * l1;
            }
        }
    } else // continuum flux is unity
        fill(&cfl[0], &cfl[n], 1.0);
}

size_t
RQ::section::selection_size() const {
    size_t j = 0;

    for (size_t i = 0; i < n; ++i)
        if (msk[i])
            ++j;

    return j;
}

double
RQ::section::cost() const {
    double a = 0.0;

    for (size_t i = 0; i < n; ++i)
        if (msk[i])
            a += res[i] * res[i];

    return 0.5 * a;
}

void
RQ::section::mask(double a, double b) {
    for (size_t i = 0; i < n; ++i)
        if (a <= wav[i] and wav[i] <= b)
            msk[i] = 0;
}

void
RQ::section::unmask(double a, double b) {
    for (size_t i = 0; i < n; ++i)
        if (a <= wav[i] and wav[i] <= b)
            msk[i] = 1;
}

void
RQ::section::integrals(double x, double fwhm, double &p, double &q) const {
    using std::erf; // since C++11
    using std::exp;

    const double c = 1.6651092223153955127063292897904020952612;
    const double d = 3.5449077018110320545963349666822903655951;
    const double b = fwhm / c;

    x /= b;

    p = 0.5 * erf(x);
    q = -(b * exp(-x * x)) / d;
}

std::istream &
RQ::section::get(std::istream &is, double a, double b) {
    using namespace std;

    const size_t room = 20000;

    vector<bool> w;
    vector<double> x;
    vector<double> y;
    vector<double> z;

    w.reserve(room);
    x.reserve(room);
    y.reserve(room);
    z.reserve(room);

    size_t i = 0;
    string line;

    while (getline(is, line) and !line.empty()) {
        // Skip comments
        if (line[0] == '#')
            continue;

        istringstream ist(line);
        bool tw;
        double tx, ty, tz;

        if (ist >> tx >> ty) {
            if (a <= tx and tx <= b) {
                x.push_back(tx);
                y.push_back(ty);

                if (ist >> tz)
                    z.push_back(tz);
                else
                    z.push_back(1.0);

                if (ist >> tw)
                    w.push_back(tw);
                else
                    w.push_back(true);

                ++i;
            }
        } else {
            is.setstate(ios_base::badbit | ios_base::failbit);

            return is;
        }
    }

    if (i > 0) {
        wav.resize(i);
        flx.resize(i);
        err.resize(i);
        msk.resize(i);

        opt.resize(i, 0.0);
        atm.resize(i, 0.0);
        cat.resize(i, 0.0);
        cfl.resize(i, 0.0);
        tfl.resize(i, 0.0);
        fit.resize(i, 0.0);
        res.resize(i, 0.0);

        n = i;

        copy(x.begin(), x.end(), &wav[0]);
        copy(y.begin(), y.end(), &flx[0]);
        copy(z.begin(), z.end(), &err[0]);
        copy(w.begin(), w.end(), &msk[0]);

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else
        is.setstate(ios_base::failbit);

    return is;
}

std::ostream &
RQ::section::put(std::ostream &os, double a, double b) const {
    using namespace std;

    if (os) {
        const int p = 8;  // precision
        const int w = 16; // width

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::fmtflags());
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < n; ++i)
            if (a <= wav[i] and wav[i] <= b) {
                const double nfl = flx[i] / cfl[i];
                const double ner = err[i] / cfl[i];
                // normalized observed flux and uncertainty

                os << setw(w) << wav[i];
                os << setw(w) << flx[i];
                os << setw(w) << err[i];
                os << setw(3) << msk[i];
                os << setw(w) << opt[i];
                os << setw(w) << atm[i];
                os << setw(w) << cat[i];
                os << setw(w) << cfl[i];
                os << setw(w) << tfl[i];
                os << setw(w) << fit[i];
                os << setw(w) << res[i];
                os << setw(w) << nfl;
                os << setw(w) << ner;

                os << '\n';
            }

        os.flush();
        os.flags(f);
    }

    return os;
}

std::istream &
RQ::operator>>(std::istream &is, std::vector<section> &s) {
    using namespace std;

    vector<section> tmp(1);

    while (is >> tmp[0])
        tmp.push_back(tmp[0]);

    if (!is.bad() and tmp.size() > 1)
        s.assign(tmp.begin() + 1, tmp.end());

    return is;
}

std::ostream &
RQ::operator<<(std::ostream &os, const std::vector<section> &s) {
    size_t i;

    for (i = 0; i + 1 < s.size(); ++i)
        os << s[i] << '\n';
    os << s[i];

    return os;
}

// References
//
// W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery (2002)
//   Numerical Recipes in C: The Art of Scientific Computing
//   Cambridge University Press, ISBN 0-521-75033-4
