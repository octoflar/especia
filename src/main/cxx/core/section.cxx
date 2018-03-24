/// @file section.cxx
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
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>

#include "section.h"

using especia::natural;
using especia::real;
using especia::sqrt_of_ln_two;
using especia::sqrt_of_pi;

especia::Section::Section()
        : wav(),
          flx(),
          unc(),
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

especia::Section::Section(size_t n_in)
        : wav(0.0, n_in),
          flx(0.0, n_in),
          unc(0.0, n_in),
          msk(true, n_in),
          opt(0.0, n_in),
          atm(0.0, n_in),
          cat(0.0, n_in),
          cfl(0.0, n_in),
          tfl(0.0, n_in),
          fit(0.0, n_in),
          res(0.0, n_in),
          n(n_in) {
}

especia::Section::Section(size_t n_in, const real x[], const real y[], const real unc[])
        : wav(x, n_in),
          flx(y, n_in),
          unc(unc, n_in),
          msk(true, n_in),
          opt(0.0, n_in),
          atm(0.0, n_in),
          cat(0.0, n_in),
          cfl(0.0, n_in),
          tfl(0.0, n_in),
          fit(0.0, n_in),
          res(0.0, n_in),
          n(n_in) {
}

especia::Section::~Section() = default;

void especia::Section::continuum(natural m, const std::valarray<real> &cat, std::valarray<real> &cfl) const {
    using std::fill;
    using std::runtime_error;
    using std::sqrt;
    using std::valarray;

    if (m > 0) {
        valarray<real> b(0.0, m);
        valarray<real> c(0.0, m);
        valarray<valarray<real>> a(b, m);
        valarray<valarray<real>> l(valarray<real>(1.0, n), m);

        if (m > 1) {
            l[1] = 2.0 * (wav - lower_bound()) / width() - 1.0;
            // Bonnet’s recursion formula
            for (natural j = 1; j + 1 < m; ++j) {
                l[j + 1] = (real(2 * j + 1) * l[1] * l[j] - real(j) * l[j - 1]) / real(j + 1);
            }
        }

        // Optimizing the background continuum is a linear optimization problem. Here the normal
        // equations are established.
        const valarray<real> p = cat / sq(unc);
        for (natural j = 0; j < m; ++j) {
            const valarray<real> &lj = l[j];
            for (natural k = j; k < m; ++k) {
                const valarray<real> &lk = l[k];
                for (size_t i = 0; i < n; ++i) {
                    if (msk[i]) {
                        a[j][k] += cat[i] * p[i] * lj[i] * lk[i];
                    }
                }
            }
            for (size_t i = 0; i < n; ++i) {
                if (msk[i]) {
                    b[j] += flx[i] * p[i] * lj[i];
                }
            }
        }
        // The normal equations are solved by means of a Cholesky decomposition (e.g. Press et al. 2002).
        //
        // W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery (2002).
        //   Numerical Recipes in C: The Art of Scientific Computing.
        //   Cambridge University Press, ISBN 0-521-75033-4.
        for (natural i = 0; i < m; ++i) {
            for (natural j = i; j < m; ++j) {
                real s = a[i][j];

                for (natural k = 0; k < i; ++k) {
                    s -= a[i][k] * a[j][k];
                }
                if (i < j) {
                    a[j][i] = s / a[i][i];
                } else if (s > 0.0) {
                    a[i][i] = sqrt(s);
                } else {
                    // The normal equations are (numerically) singular.
                    throw runtime_error(
                            "especia::section::continuum(): Error: normal equations are numerically singular");
                }
            }
        }
        for (natural i = 0; i < m; ++i) {
            real s = b[i];

            for (natural k = 0; k < i; ++k) {
                s -= a[i][k] * c[k];
            }

            c[i] = s / a[i][i];
        }
        for (natural i = m - 1; i + 1 > 0; --i) {
            real s = c[i];

            for (natural k = i + 1; k < m; ++k) {
                s -= a[k][i] * c[k];
            }

            c[i] = s / a[i][i];
        }

        // Compute the continuum flux. The first Legendre term is a constant.
        cfl = c[0];
        // The other terms depend on the abcissa value.
        for (natural k = 1; k < m; ++k) {
            cfl += c[k] * l[k];
        }
    } else {
        fill(&cfl[0], &cfl[n], 1.0);
    }
}

size_t especia::Section::valid_data_count() const {
    size_t count = 0;

    for (size_t i = 0; i < n; ++i) {
        if (msk[i]) {
            ++count;
        }
    }

    return count;
}

real especia::Section::cost() const {
    real cost = 0.0;

    for (size_t i = 0; i < n; ++i) {
        if (msk[i]) {
            cost += res[i] * res[i];
        }
    }

    return 0.5 * cost;
}

void especia::Section::mask(real a, real b) {
    for (size_t i = 0; i < n; ++i) {
        if (a <= wav[i] and wav[i] <= b) {
            msk[i] = false;
        }
    }
}

void especia::Section::primitive(const real &x, const real &h, real &p, real &q) {
    using std::erf; // C++11
    using std::exp;

    const real b = h / sqrt_of_ln_two;
    const real d = b / sqrt_of_pi;

    p = 0.5 * erf(x / b);
    q = 0.5 * exp(-sq(x / b)) * (-d);
}

void especia::Section::supersample(const std::valarray<real> &source, natural k, std::valarray<real> &target) {
    for (natural is = 0, it = 0; is < source.size(); ++is, it += k) {
        target[it] = source[is];
    }
    for (natural j = 1; j < k; ++j) {
        const real w = real(j) / real(k);

        for (size_t is = 0, it = j; is + 1 < source.size(); ++is, it += k) {
            target[it] = source[is] + w * (source[is + 1] - source[is]);
        }
    }
}

std::istream &especia::Section::get(std::istream &is, real a, real b) {
    using namespace std;

    const size_t room = 20000;

    vector<bool> w;
    vector<real> x;
    vector<real> y;
    vector<real> z;

    w.reserve(room);
    x.reserve(room);
    y.reserve(room);
    z.reserve(room);

    size_t i = 0;
    string line;

    while (getline(is, line) and !line.empty()) {
        // Skip comments.
        if (line[0] == '#' or line[0] == '%' or line[0] == '!') {
            continue;
        }

        istringstream ist(line);
        bool tw;
        real tx, ty, tz;

        if (ist >> tx >> ty) {
            if (a <= tx and tx <= b) {
                x.push_back(tx);
                y.push_back(ty);

                if (ist >> tz) {
                    z.push_back(tz);
                } else {
                    z.push_back(1.0);
                }
                if (ist >> tw) {
                    w.push_back(tw);
                } else {
                    w.push_back(true);
                }

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
        unc.resize(i);
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
        copy(z.begin(), z.end(), &unc[0]);
        copy(w.begin(), w.end(), &msk[0]);

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else {
        is.setstate(ios_base::failbit);
    }

    return is;
}

std::ostream &especia::Section::put(std::ostream &os, real a, real b) const {
    using namespace std;

    if (os) {
        // The precision.
        const natural p = 8;
        // The width of the output field.
        const natural w = 16;

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::fmtflags());
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < n; ++i)
            if (a <= wav[i] and wav[i] <= b) {
                // The normalized observed spectral flux and its uncertainty.
                const real nfl = flx[i] / cfl[i];
                const real nun = unc[i] / cfl[i];

                os << setw(w) << wav[i]; // 1
                os << setw(w) << flx[i]; // 2
                os << setw(w) << unc[i]; // 3
                os << setw(3) << msk[i]; // 4
                os << setw(w) << opt[i]; // 5
                os << setw(w) << atm[i]; // 6
                os << setw(w) << cat[i]; // 7
                os << setw(w) << cfl[i]; // 8
                os << setw(w) << tfl[i]; // 9
                os << setw(w) << fit[i]; // 10
                os << setw(w) << res[i]; // 11
                os << setw(w) << nfl;    // 12
                os << setw(w) << nun;    // 13

                os << '\n';
            }

        os.flush();
        os.flags(f);
    }

    return os;
}

std::istream &especia::operator>>(std::istream &is, std::vector<Section> &sections) {
    using std::vector;

    Section s;

    while (is >> s) {
        sections.push_back(s);
    }

    return is;
}

std::ostream &especia::operator<<(std::ostream &os, const std::vector<Section> &sections) {
    size_t i;

    for (i = 0; i + 1 < sections.size(); ++i) {
        os << sections[i] << '\n';
    }
    os << sections[i];

    return os;
}
