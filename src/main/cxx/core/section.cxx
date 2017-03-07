/// @file section.cxx
/// Class for modeling spectroscopic data sections.
/// Copyright (c) 2016 Ralf Quast
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

using especia::N_elem;
using especia::R_elem;

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

especia::Section::Section(size_t m)
        : wav(0.0, m),
          flx(0.0, m),
          unc(0.0, m),
          msk(true, m),
          opt(0.0, m),
          atm(0.0, m),
          cat(0.0, m),
          cfl(0.0, m),
          tfl(0.0, m),
          fit(0.0, m),
          res(0.0, m),
          n(m) {
}

especia::Section::Section(size_t m, const R_elem x[], const R_elem y[], const R_elem unc[])
        : wav(x, m),
          flx(y, m),
          unc(unc, m),
          msk(true, m),
          opt(0.0, m),
          atm(0.0, m),
          cat(0.0, m),
          cfl(0.0, m),
          tfl(0.0, m),
          fit(0.0, m),
          res(0.0, m),
          n(m) {
}

especia::Section::~Section() {
}

void especia::Section::continuum(N_elem m, const R_elem cat[], R_elem cfl[]) const throw(std::runtime_error) {
    using std::fill;
    using std::runtime_error;
    using std::sqrt;
    using std::valarray;

    if (m > 0) {
        valarray<R_elem> b(0.0, m);
        valarray<R_elem> c(0.0, m);
        valarray<R_elem> l(1.0, m);

        valarray<valarray<R_elem> > a(b, m);

        // Optimizing the background continuum is a linear optimization problem. Here the normal
        // equations are established.
        for (size_t i = 0; i < n; ++i) {
            if (msk[i]) {
                // Map the wavelengths onto the interval [-1, 1]
                const R_elem x = 2.0 * (wav[i] - wav[0]) / width() - 1.0;

                R_elem l1 = 1.0;
                R_elem l2 = 0.0;

                // Compute the higher-order Legendre basis polynomials.
                for (N_elem j = 1; j < m; ++j) {
                    const R_elem l3 = l2;

                    l2 = l1;
                    l1 = ((2 * j - 1) * x * l2 - (j - 1) * l3) / j;
                    l[j] = l1;
                }
                // Establish the normal equations.
                for (N_elem j = 0; j < m; ++j) {
                    for (N_elem k = j; k < m; ++k) {
                        a[j][k] += (cat[i] * cat[i] * l[j] * l[k]) / (unc[i] * unc[i]);
                    }

                    b[j] += (flx[i] * cat[i] * l[j]) / (unc[i] * unc[i]);
                }
            }
        }
        // The normal equations are solved by means of a Cholesky decomposition (e.g. Press et al. 2002).
        //
        // W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery (2002).
        //   Numerical Recipes in C: The Art of Scientific Computing.
        //   Cambridge University Press, ISBN 0-521-75033-4.
        for (N_elem i = 0; i < m; ++i) {
            for (N_elem j = i; j < m; ++j) {
                R_elem s = a[i][j];

                for (N_elem k = 0; k < i; ++k) {
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
        for (N_elem i = 0; i < m; ++i) {
            R_elem s = b[i];

            for (N_elem k = 0; k < i; ++k) {
                s -= a[i][k] * c[k];
            }

            c[i] = s / a[i][i];
        }
        for (N_elem i = m - 1; i + 1 > 0; --i) {
            R_elem s = c[i];

            for (N_elem k = i + 1; k < m; ++k) {
                s -= a[k][i] * c[k];
            }

            c[i] = s / a[i][i];
        }

        // Compute the continuum flux.
        for (size_t i = 0; i < n; ++i) {
            const R_elem x = 2.0 * (wav[i] - wav[0]) / width() - 1.0;

            R_elem l1 = 1.0;
            R_elem l2 = 0.0;

            cfl[i] = c[0];

            for (N_elem k = 1; k < m; ++k) {
                const R_elem l3 = l2;

                l2 = l1;
                l1 = ((2 * k - 1) * x * l2 - (k - 1) * l3) / k;
                cfl[i] += c[k] * l1;
            }
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

R_elem especia::Section::cost() const {
    R_elem cost = 0.0;

    for (size_t i = 0; i < n; ++i) {
        if (msk[i]) {
            cost += res[i] * res[i];
        }
    }

    return 0.5 * cost;
}

void especia::Section::mask(R_elem a, R_elem b) {
    for (size_t i = 0; i < n; ++i) {
        if (a <= wav[i] and wav[i] <= b) {
            msk[i] = 0;
        }
    }
}

void especia::Section::primitive(R_elem x, R_elem h, R_elem &p, R_elem &q) const {
    using std::erf; // C++11
    using std::exp;

    const R_elem c = 1.6651092223153955127063292897904020952612;
    const R_elem d = 3.5449077018110320545963349666822903655951;
    const R_elem b = 2.0 * h / c;

    x /= b;

    p = 0.5 * erf(x);
    q = -(b * exp(-x * x)) / d;
}

std::istream &especia::Section::get(std::istream &is, R_elem a, R_elem b) {
    using namespace std;

    const size_t room = 20000;

    vector<bool> w;
    vector<R_elem> x;
    vector<R_elem> y;
    vector<R_elem> z;

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
        R_elem tx, ty, tz;

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

std::ostream &especia::Section::put(std::ostream &os, R_elem a, R_elem b) const {
    using namespace std;

    if (os) {
        // The precision.
        const N_elem p = 8;
        // The width of the output field.
        const N_elem w = 16;

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::fmtflags());
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < n; ++i)
            if (a <= wav[i] and wav[i] <= b) {
                // The normalized observed spectral flux and its uncertainty.
                const R_elem nfl = flx[i] / cfl[i];
                const R_elem nun = unc[i] / cfl[i];

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
