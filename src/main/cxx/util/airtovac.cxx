// Utility: convert wavelength in spectroscopic data from air to vacuum
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
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

#include "../base.h"

using namespace std;

istream &get(istream &is, valarray<double> &x, valarray<double> &y, valarray<double> &z, int skip = 0) {
    const size_t room = 20000;

    vector<double> u;
    vector<double> v;
    vector<double> w;

    u.reserve(room);
    v.reserve(room);
    w.reserve(room);

    size_t n = 0;
    string s;

    while (getline(is, s))
        if (skip <= 0) {
            istringstream ist(s);
            double a, b, c;

            if (ist >> a >> b) {
                u.push_back(a);
                v.push_back(b);
                if (ist >> c)
                    w.push_back(c);

                ++n;
            } else {
                is.setstate(ios_base::badbit | ios_base::failbit);

                return is;
            }
        } else
            --skip;

    if (n > 0 and is.eof()) {
        x.resize(n);
        y.resize(n);
        z.resize(n);

        copy(u.begin(), u.end(), &x[0]);
        copy(v.begin(), v.end(), &y[0]);
        copy(w.begin(), w.end(), &z[0]);

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else
        is.setstate(ios_base::failbit);

    return is;
}

ostream &put(ostream &os, const valarray<double> &x, const valarray<double> &y, const valarray<double> &z) {
    if (os) {
        const int p = 6;  // precision
        const int w = 14; // width

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < x.size(); ++i) {
            os.setf(ios_base::fixed, ios_base::floatfield);
            os << setw(w) << x[i];
            os.setf(ios_base::scientific, ios_base::floatfield);
            os << setw(w) << y[i];
            if (z.size() > 0)
                os << setw(w) << z[i];
            os << '\n';
        }

        os.flush();
        os.flags(f);
    }

    return os;
}

void vactoair(double x, double &y, double &z) {
    using especia::sqr;

/*  const double a = 1.0000643280 + 2.5540e-10 / (0.0000410 - x * x) + 2.949810e-08 / (0.000146 - x * x);
        // Edlen (1953) */
    const double a = 1.0000834213 + 1.5997e-10 / (0.0000389 - x * x) + 2.406030e-08 / (0.000130 - x * x);
        // Edlen (1966)

    y = a * x;
/*  z = a + x * ((5.1080e-10 * x) / sqr(0.0000410 - x * x) + (5.89962e-08 * x) / sqr(0.000146 - x * x));
        // Edlen (1953), the first derivative of y with respect to x */
    z = a + x * ((3.1994e-10 * x) / sqr(0.0000389 - x * x) + (4.81206e-08 * x) / sqr(0.000130 - x * x));
        // Edlen (1966), the first derivative of y with respect to x
}

double find_root(void f(double, double &, double &), double c, double x, double accuracy_goal) throw(runtime_error) {
    double d, y, z;
    unsigned i = 0;

    do {
        f(x, y, z);
        d = (y - c) / z;
        x -= d;
    } while (++i < 100 and abs(d) > accuracy_goal * x);

    if (i == 100)
        throw runtime_error("find_root(): Error: iteration stopped");

    return x;
}

/**
 * Utility to convert wavelength in spectroscopic data from air
 * to vacuum.
 *
 * @remark Usage: airtovac [SKIP] [ACCURACY] < ISTREAM > OSTREAM
 *
 * Further reading:
 *
 * B. Edlen (1966).
 *   *The refractive index of air.*
 *   Metrologia, 2, 2, 71-80.
 *   http://dx.doi.org/10.1088/0026-1394/2/2/002
 *
 * B. Edlen (1953).
 *   *The dispersion of standard air.*
 *   Journal of the Optical Society of America, 43, 5, 339.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The number of lines to skip (optional, default = 0).
 *
 * @c argv[2] The accuracy goal (optional, dafault = 1.0E-8).
 * @endparblock
 * @return an exit code.
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    int skip = 0;
    double accuracy_goal = 1.0e-08;

    if (argc == 3) {
        accuracy_goal = atof(argv[2]);
    }
    if (argc == 3 || argc == 2) {
        skip = atoi(argv[1]);
    }
    if (argc == 3 || argc == 2 || argc == 1) {
        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        if (get(cin, x, y, z, skip)) {
            try {
                for (size_t i = 0; i < x.size(); ++i)
                    x[i] = 10.0 / find_root(vactoair, 10.0 / x[i], 10.0 / x[i], accuracy_goal);
            }
            catch (exception &e) {
                cerr << pname << ": " << e.what() << endl;
                return 3;
            }

            put(cout, x, y, z);
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " [SKIP] [ACCURACY] < ISTREAM > OSTREAM" << endl;
        return 1;
    }
}
