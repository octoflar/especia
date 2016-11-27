// Utility: convert wavelength in spectroscopic data from vacuum to air
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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

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
            os << setw(w) << z[i];
            os << '\n';
        }

        os.flush();
        os.flags(f);
    }

    return os;
}

double
vactoair(double x) {
/*  const double a = 1.0000643280 + 2.5540e-10 / (0.0000410 - x * x) + 2.949810e-08 / (0.000146 - x * x);
        // Edlen (1953) */
    const double a = 1.0000834213 + 1.5997e-10 / (0.0000389 - x * x) + 2.406030e-08 / (0.000130 - x * x);
        // Edlen (1966)

    return a * x;
}

int main(int argc, char *argv[]) {
    const char *pname = argv[0];
    int skip = 0;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "missing_default_case"
    switch (argc) {
        case 2:
            skip = atoi(argv[1]);
        case 1:
            valarray<double> x;
            valarray<double> y;
            valarray<double> z;

            if (get(cin, x, y, z, skip)) {
                for (size_t i = 0; i < x.size(); ++i)
                    x[i] = 10.0 / vactoair(10.0 / x[i]);

                put(cout, x, y, z);
            } else {
                cerr << pname << ": input failure" << endl;
                return 2;
            }

            return 0;
        default:
            cout << "usage: " << pname << " [SKIP] < ISTREAM > OSTREAM" << endl;
            return 1;
    }
}

// References
//
// B. Edlen (1966)
//   The refractive index of air
//   Metrologia, 2, 2, 71-80
//   http://dx.doi.org/10.1088/0026-1394/2/2/002
//
// B. Edlen (1953)
//   The dispersion of standard air
//   Journal of the Optical Society of America, 43, 5, 339
