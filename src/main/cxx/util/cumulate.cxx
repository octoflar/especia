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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>

#include "../base.h"

using namespace std;


class Frame {
public:
    friend class Frame_Stack;

    Frame() : x(), y(), z(), c(), d(), i(0), j(0), n(0) {
    }

    ~Frame() {
    }

    Frame &resize(size_t m) {
        x.resize(m, 0.0);
        y.resize(m, 0.0);
        z.resize(m, 0.0);
        c.resize(m, 0.0);
        d.resize(m, 0.0);

        i = 0;
        j = m - 1;
        n = m;

        return *this;
    }

    Frame &resample(double r = 0.0) {
        using especia::kilo;

        const double h = (r > 0.0) ? 0.5 * center() / (r * kilo) : width() / (n - 1);
        const size_t m = (r > 0.0) ? static_cast<size_t>(width() / h) + 1 : n;

        valarray<double> u(m);
        valarray<double> v(m);
        valarray<double> w(m);

        for (size_t i = 0; i < m; ++i) {
            splint(u[i] = lower_bound() + i * h, v[i], w[i]);
        }
        resize(m);

        copy(&u[0], &u[n], &x[0]);
        copy(&v[0], &v[n], &y[0]);
        copy(&w[0], &w[n], &z[0]);

        spline();

        return *this;
    }

    Frame &scale(double a) {
        y *= a;
        z *= a;
        spline();

        return *this;
    }

    double lower_bound() const {
        return x[0];
    }

    double upper_bound() const {
        return (n > 1) ? x[n - 1] : x[0];
    }

    double center() const {
        return 0.5 * (lower_bound() + upper_bound());
    }

    double width() const {
        return upper_bound() - lower_bound();
    }

    double median() const {
        using especia::Indirect_Compare;

        valarray<size_t> index(n);
        for (size_t i = 0; i < n; ++i)
            index[i] = i;

        nth_element(&index[n / 3], &index[n >> 1], &index[n / 3 << 1],
                    Indirect_Compare<double, less<double> >(y, less<double>()));
        cout << "frame::median(): Message: median is " << y[index[n >> 1]] << endl;

        return y[index[n >> 1]];
    }

    size_t data_count() const {
        return n;
    }

    istream &get(istream &is, double a = 0.0, double b = numeric_limits<double>::max()) {
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

    ostream &put(ostream &os, double a = 0.0, double b = numeric_limits<double>::max()) const {
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

private:
    valarray<double> x; // wavelength data
    valarray<double> y; // flux data
    valarray<double> z; // flux data uncertainty
    valarray<double> c; // flux cubic spline data
    valarray<double> d; // uncertainty cubic spline data

    mutable size_t i;
    mutable size_t j;
    size_t n; // number of data points

    void spline() {
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

    void splint(double u, double &v, double &w) const throw(runtime_error) {
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
        } else {
            throw runtime_error("frame::splint(): Error: bad abscissa table");
        }
    }
};


class Frame_Stack {
public:
    Frame_Stack() : frames() {
    }

    ~Frame_Stack() {
    }

    Frame_Stack &align() {
        return *this;
    }

    Frame combine() const {
        Frame f;

        if (!frames.empty()) {
            f.resize(frames.front().data_count());

            for (size_t i = 0; i < f.data_count(); ++i) {
                f.x[i] = frames.front().x[i];
                f.y[i] = 0.0;
                f.z[i] = 0.0;

                for (vector<Frame>::const_iterator j = frames.begin(); j < frames.end(); ++j) {
                    if (j->lower_bound() <= f.x[i] and f.x[i] <= j->upper_bound()) {
                        double w, y, z;

                        j->splint(f.x[i], y, z);

                        if (z > 0.0) {
                            w = 1.0 / (z * z);
                            f.y[i] += w * y;
                            f.z[i] += w;
                        }
                    }
                }
                if (f.z[i] > 0.0) {
                    f.y[i] /= f.z[i];
                    f.z[i] = sqrt(1.0 / f.z[i]);
                }
            }
            f.spline();
        }
        return f;
    }

    Frame_Stack &scale() {
        using especia::Indirect_Compare;

        valarray<double> m(size());
        valarray<size_t> j(size());

        for (size_t i = 0; i < size(); ++i) {
            m[i] = frames[i].median();
            j[i] = i;
        }
        nth_element(&j[0], &j[size() >> 1], &j[size()], Indirect_Compare<double, less<double> >(m, less<double>()));

        const size_t k = j[size() >> 1];

        for (size_t i = 0; i < size(); ++i) {
            if (i != k) {
                const double a = m[k] / m[i];

                frames[i].scale(a);
                clog << "stack::rescale(): Message: frame rescaled by " << a << endl;
            }
        }
        return *this;
    }

    size_t size() const {
        return frames.size();
    }

    vector<Frame> frames;
};


istream &operator>>(istream &is, Frame &frame) {
    return frame.get(is);
}

ostream &operator<<(ostream &os, const Frame &frame) {
    return frame.put(os);
}

istream &operator>>(istream &is, Frame_Stack &stack) {
    using namespace std;

    Frame frame;

    while (is >> frame)
        stack.frames.push_back(frame);

    if (!is.bad() and stack.size() > 0) {
        is.clear(is.rdstate() & ~ios_base::failbit);
    }

    return is;
}

ostream &operator<<(ostream &os, const Frame_Stack &stack) {
    size_t i;

    for (i = 0; i + 1 < stack.size(); ++i)
        os << stack.frames[i] << '\n';
    os << stack.frames[i];

    return os;
}

/**
 * Utility to accumulate multiple spectra of the same target,
 * acquired at different dates.
 *
 * The utility (1) reads a stack of spectroscopic data frames
 * from standard input, where the frame separator is an empty
 * line, (2) scales all frames are to the same median flux level,
 * (3) computes the weighted average of all frames, and (4) samples
 * the weighted average frame to equidistant wavelengths.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The spectral resolution of the data (1.0E+3). This parameter is optional.
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: cumulate [RESOLUTION] < ISTREAM > OSTREAM
 *
 * @deprecated The current implementation interprets the spectroscopic data as a cubic
 * spline, therefore the frames need not be co-aligned. The final resampling step might be
 * sufficiently accurate, but does not necessarily conserve flux. A flux-conserving
 * co-alignment of frames before averaging will be more accurate than the procedure
 * implemented here, the use of which is discouraged. My recommended way to combine
 * spectroscopic data is to define a model definition block for each individual exposure
 * taken.
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    double resolution = 0.0;

    if (argc == 2) {
        resolution = atof(argv[1]);
    }
    if (argc == 2 or argc == 1) {
        Frame_Stack stack;

        if (cin >> stack) {
            try {
                cout << stack.scale().align().combine().resample(resolution);
            }
            catch (exception &e) {
                cerr << pname << ": " << e.what() << endl;
                return 3;
            }
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " [RESOLUTION] < ISTREAM > OSTREAM" << endl;
        return 1;
    }
}
