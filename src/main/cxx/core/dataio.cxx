/// @file dataio.cxx
/// Data input and output procedures.
/// Copyright (c) 2021 Ralf Quast
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
#include <iomanip>
#include <sstream>
#include <vector>

#include "dataio.h"

using namespace std;

istream &especia::get(istream &is, valarray<real> &x, valarray<real> &y, natural skip) {
    const size_t room = 20000;

    vector<real> u;
    vector<real> v;

    u.reserve(room);
    v.reserve(room);

    size_t n = 0;
    string s;

    while (getline(is, s)) {
        if (skip <= 0) {
            istringstream ist(s);
            real a, b;

            if (ist >> a >> b) {
                u.push_back(a);
                v.push_back(b);

                ++n;
            } else {
                is.setstate(ios_base::badbit | ios_base::failbit);

                return is;
            }
        } else {
            --skip;
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

istream &especia::get(istream &is, valarray<real> &x, valarray<real> &y, valarray<real> &z, natural skip) {
    const size_t room = 20000;

    vector<real> u;
    vector<real> v;
    vector<real> w;

    u.reserve(room);
    v.reserve(room);
    w.reserve(room);

    size_t n = 0;
    string s;

    while (getline(is, s)) {
        if (skip <= 0) {
            istringstream ist(s);
            real a, b, c;

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
        } else {
            --skip;
        }
    }

    if (n > 0 and is.eof()) {
        x.resize(n);
        y.resize(n);
        z.resize(n);

        copy(u.begin(), u.end(), &x[0]);
        copy(v.begin(), v.end(), &y[0]);
        copy(w.begin(), w.end(), &z[0]);

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else {
        is.setstate(ios_base::failbit);
    }

    return is;
}

ostream &especia::put(ostream &os, const valarray<real> &x, const valarray<real> &y, const valarray<real> &z) {
    if (os) {
        const natural p = 6;  // precision
        const natural w = 14; // width

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < x.size(); ++i) {
            os.setf(ios_base::fixed, ios_base::floatfield);
            os << setw(w) << x[i];
            os.setf(ios_base::scientific, ios_base::floatfield);
            os << setw(w) << y[i];
            if (z.size() > 0) {
                os << setw(w) << z[i];
            }
            os << '\n';
        }

        os.flush();
        os.flags(f);
    }

    return os;
}
