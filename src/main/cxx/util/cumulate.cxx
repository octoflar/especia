// Utility: accumulate multiple spectra (weighted average)
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
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

using namespace std;

template<class number, class comparation>
class indirect_comparation;

class frame;

class stack;

template<class number, class comparation>
class indirect_comparation {
public:
    indirect_comparation(const valarray<number> &y, const comparation &c): x(y), comp(c) {
    }

    ~indirect_comparation() {
    }

    // comparation operator
    bool operator()(size_t i, size_t j) {
        return comp(x[i], x[j]);
    }

private:
    const valarray<number> &x;
    const comparation &comp;
};


class frame {
public:
    friend class stack;

    frame();

    frame(size_t n);

    frame(const double x[], const double y[], const double z[], size_t n);

    ~frame();

    void resize(size_t m);

    void resample(double r = 0.0);

    void scale(double a);

    void operator()(double x, double &y, double &z) const;

    double blu_end() const;

    double red_end() const;

    double center() const;

    double length() const;

    double median() const;

    size_t size() const;

    istream &get(istream &is,
                 double a = 0.0, double b = numeric_limits<double>::max());

    ostream &put(ostream &os,
                 double a = 0.0, double b = numeric_limits<double>::max()) const;

private:
    valarray<double> x; // wavelength data
    valarray<double> y; // flux data
    valarray<double> z; // flux data uncertainty
    valarray<double> c; // flux cubic spline data
    valarray<double> d; // uncertainty cubic spline data

    mutable size_t i;
    mutable size_t j;
    size_t n; // number of data points

    void spline();

    void splint(double x, double &y, double &z) const throw(runtime_error);
};

frame::frame()
        : x(),
          y(),
          z(),
          c(),
          d(),
          i(0),
          j(0),
          n(0) {
}

frame::frame(size_t m)
        : x(0.0, m),
          y(0.0, m),
          z(0.0, m),
          c(0.0, m),
          d(0.0, m),
          i(0),
          j(m - 1),
          n(m) {
}

frame::frame(const double u[], const double v[], const double w[], size_t m)
        : x(u, m),
          y(v, m),
          z(w, m),
          c(0.0, m),
          d(0.0, m),
          i(0),
          j(m - 1),
          n(m) {
    spline();
}

frame::~frame() {
}

void frame::resize(size_t m) {
    x.resize(m, 0.0);
    y.resize(m, 0.0);
    z.resize(m, 0.0);
    c.resize(m, 0.0);
    d.resize(m, 0.0);

    i = 0;
    j = m - 1;
    n = m;
}

void frame::resample(double r) {
    const double h = (r > 0.0) ? center() / (2.0 * r) : length() / (n - 1);
    const size_t m = (r > 0.0) ? static_cast<size_t>(length() / h) + 1 : n;

    valarray<double> u(m);
    valarray<double> v(m);
    valarray<double> w(m);

    for (size_t i = 0; i < m; ++i)
        splint(u[i] = blu_end() + i * h, v[i], w[i]);
    resize(n = m);

    copy(&u[0], &u[n], &x[0]);
    copy(&v[0], &v[n], &y[0]);
    copy(&w[0], &w[n], &z[0]);

    spline();
}

void frame::scale(double a) {
    y *= a;
    z *= a;
    spline();
}

inline
void frame::operator()(double x, double &y, double &z) const {
    splint(x, y, z);
}

inline
double frame::blu_end() const {
    return x[0];
}

inline
double frame::red_end() const {
    return (n > 1) ? x[n - 1] : x[0];
}

inline
double frame::center() const {
    return 0.5 * (blu_end() + red_end());
}

inline
double frame::length() const {
    return red_end() - blu_end();
}

double frame::median() const {
    valarray<size_t> index(n);
    for (size_t i = 0; i < n; ++i)
        index[i] = i;

    nth_element(&index[n / 3], &index[n >> 1], &index[n / 3 << 1], indirect_comparation<double,
            less<double> >(y, less<double>()));
    clog << "frame::median(): Message: median is " << y[index[n >> 1]] << endl;

    return y[index[n >> 1]];
}

inline
size_t frame::size() const {
    return n;
}

istream &frame::get(istream &is, double a, double b) {
    const size_t room = 20000;

    vector<double> u;
    vector<double> v;
    vector<double> w;

    u.reserve(room);
    v.reserve(room);
    w.reserve(room);

    size_t m = 0;
    string line;

    while (getline(is, line) and !line.empty()) {
        istringstream ist(line);
        double x, y, z;

        if (ist >> x >> y >> z) {
            if (a <= x and x <= b) {
                u.push_back(x);
                v.push_back(y);
                w.push_back(z);

                ++m;
            }
        } else {
            is.setstate(ios_base::badbit | ios_base::failbit);

            return is;
        }
    }

    if (m > 0) {
        x.resize(m);
        y.resize(m);
        z.resize(m);
        c.resize(m);
        d.resize(m);

        i = 0;
        j = m - 1;
        n = m;

        copy(u.begin(), u.end(), &x[0]);
        copy(v.begin(), v.end(), &y[0]);
        copy(w.begin(), w.end(), &z[0]);

        spline();

        is.clear(is.rdstate() & ~ios_base::failbit);
    } else
        is.setstate(ios_base::failbit);

    return is;
}

ostream &frame::put(ostream &os, double a, double b) const {
    if (os) {
        const int p = 6;  // precision
        const int w = 14; // width

        const ios_base::fmtflags f = os.flags();

        os.setf(ios_base::fmtflags());
        os.setf(ios_base::right, ios_base::adjustfield);
        os.precision(p);

        for (size_t i = 0; i < n; ++i)
            if (a <= x[i] and x[i] <= b) {
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

void frame::spline() {
    valarray<double> u(0.0, n);
    valarray<double> v(0.0, n);

    for (size_t i = 1; i + 1 < n; ++i) {
        const double s = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        const double p = s * c[i - 1] + 2.0;
        const double q = s * d[i - 1] + 2.0;

        c[i] = (s - 1.0) / p;
        d[i] = (s - 1.0) / q;
        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - s * u[i - 1]) / p;
        v[i] = (z[i + 1] - z[i]) / (x[i + 1] - x[i]) - (z[i] - z[i - 1]) / (x[i] - x[i - 1]);
        v[i] = (6.0 * v[i] / (x[i + 1] - x[i - 1]) - s * v[i - 1]) / q;
    }

    for (size_t i = n - 2; i > 1; --i) {
        c[i] = c[i] * c[i + 1] + u[i];
        d[i] = d[i] * d[i + 1] + v[i];
    }
}

void frame::splint(double u, double &v, double &w) const throw(runtime_error) {
    if (i > 0 and x[--i] > u)
        i = 0;
    if (j + 1 < n and x[++j] < u)
        j = n - 1;

    while (j > i + 1) {
        const size_t k = (i + j) >> 1;

        if (x[k] > u)
            j = k;
        else
            i = k;
    }

    if (const double h = x[j] - x[i]) {
        const double a = (x[j] - u) / h;
        const double b = (u - x[i]) / h;

        v = a * y[i] + b * y[j];
        w = a * z[i] + b * z[j];
        // linear interpolation
        v += ((a * a * a - a) * y[i] + (b * b * b - b) * y[j]) * (h * h) / 6.0;
        w += ((a * a * a - a) * z[i] + (b * b * b - b) * z[j]) * (h * h) / 6.0;
        // cubic spline interpolation
    } else
        throw runtime_error("frame::splint(): Error: bad abscissa table");
}

class stack {
public:
    stack();

    ~stack();

    void coadd(frame &f) const;

    void align();

    void rescale();

    size_t size() const;

    vector<frame> frames;
};

stack::stack()
        : frames() {
}

stack::~stack() {
}

void stack::coadd(frame &f) const {
    if (!frames.empty()) {
        f.resize(frames.front().size());

        for (size_t i = 0; i < f.size(); ++i) {
            f.x[i] = frames.front().x[i];
            f.y[i] = 0.0;
            f.z[i] = 0.0;

            for (vector<frame>::const_iterator j = frames.begin(); j < frames.end(); ++j)
                if (j->blu_end() <= f.x[i] and f.x[i] <= j->red_end()) {
                    double w, y, z;

                    j->splint(f.x[i], y, z);

                    if (z > 0.0) {
                        w = 1.0 / (z * z);
                        f.y[i] += w * y;
                        f.z[i] += w;
                    }
                }

            if (f.z[i] > 0.0) {
                f.y[i] /= f.z[i];
                f.z[i] = sqrt(1.0 / f.z[i]);
            }
        }
        f.spline();
    }
}

void stack::align() {
    return; // intentionally do nothing
}

void stack::rescale() {
    valarray<double> m(size());
    valarray<size_t> j(size());

    for (size_t i = 0; i < size(); ++i) {
        m[i] = frames[i].median();
        j[i] = i;
    }
    nth_element(&j[0], &j[size() >> 1], &j[size()], indirect_comparation<double, less<double> >(m, less<double>()));

    const size_t k = j[size() >> 1];

    for (size_t i = 0; i < size(); ++i)
        if (i != k) {
            const double a = m[k] / m[i];

            frames[i].scale(a);
            clog << "stack::rescale(): Message: frame rescaled by " << a << endl;
        }
}

inline
size_t stack::size() const {
    return frames.size();
}

inline
istream &operator>>(istream &is, frame &f) {
    return f.get(is);
}

inline
ostream &operator<<(ostream &os, const frame &f) {
    return f.put(os);
}

istream &operator>>(istream &is, stack &s) {
    using namespace std;

    frame f;

    while (is >> f)
        s.frames.push_back(f);

    if (!is.bad() and s.size() > 1)
        is.clear(is.rdstate() & ~ios_base::failbit);

    return is;
}

ostream &operator<<(ostream &os, const stack &s) {
    size_t i;

    for (i = 0; i + 1 < s.size(); ++i)
        os << s.frames[i] << '\n';
    os << s.frames[i];

    return os;
}

//
// The basic procedure is:
//
// 1. A stack of spectroscopic data frames is read from standard input. The frame
//    separator is an empty line.
// 2. All frames are rescaled to the same median flux level.
// 3. The weighted average of all frames is computed.
// 4. The weighted average frame is resampled to equidistant wavelengths.
//
// The current implementation interprets the spectroscopic data as a cubic spline,
// therefore the frames need not be co-aligned. The final resampling step might be
// sufficiently accurate, but does not necessarily conserve flux.
//
// A flux-conserving co-alignment of frames before averaging might be more accurate
// than the procedure implemented here.
//
int main(int argc, char *argv[]) {
    const char *pname = argv[0];
    double resolution = 0.0;

    switch (argc) {
        case 2:
            resolution = atof(argv[1]);
        case 1:
            stack s;

            if (cin >> s) {
                frame f;

                try {
                    s.rescale();
                    s.align();
                    s.coadd(f);
                    f.resample(resolution);
                }
                catch (exception &e) {
                    cerr << pname << ": " << e.what() << endl;
                    return 0;
                }

                cout << f;
            } else {
                cerr << pname << ": input failure" << endl;
                return 2;
            }

            return 0;
        default:
            cout << "usage: " << pname << " [RESOLUTION] < ISTREAM > OSTREAM" << endl;
            return 1;
    }
}
