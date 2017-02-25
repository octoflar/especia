// Utility: merge separate spectral flux and uncertainty data to three-column format
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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

using namespace std;

istream &get(istream &is, valarray<double> &x, valarray<double> &y) {
    const size_t room = 20000;

    vector<double> u;
    vector<double> v;

    u.reserve(room);
    v.reserve(room);

    size_t n = 0;
    string s;

    while (getline(is, s)) {
        istringstream ist(s);
        double a, b;

        if (ist >> a >> b) {
            u.push_back(a);
            v.push_back(b);

            ++n;
        } else {
            is.setstate(ios_base::badbit | ios_base::failbit);

            return is;
        }
    }

    if (n > 0 and is.eof()) {
        x.resize(n);
        y.resize(n);

        copy(u.begin(), u.end(), &x[0]);
        copy(v.begin(), v.end(), &y[0]);

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else {
        is.setstate(ios_base::failbit);
    }

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

/**
 * Utility to merge separate spectral flux and uncertainty data to
 * three-column format.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The path name of the flux data file.
 *
 * @c argv[2] The path name of the flux uncertainty data file.
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: threecol FLUX UNCERTAINTY > OSTREAM
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    if (argc == 3) {
        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        ifstream fxy(argv[1]);
        ifstream fxz(argv[2]);

        get(fxy, x, y);
        get(fxz, x, z);

        fxy.close();
        fxz.close();

        if (fxy and fxz and y.size() == z.size())
            put(cout, x, y, z);
        else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " FLUX UNCERTAINTY > OSTREAM" << endl;
        return 1;
    }
}
