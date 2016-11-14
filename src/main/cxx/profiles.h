// Profile funtions
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
#ifndef RQ_PROFILE_FUNCTIONS_H
#define RQ_PROFILE_FUNCTIONS_H

#include <cmath>
#include <cstddef>
#include <vector>

namespace RQ {
    // Mathematical constants
    const double ln_of_2 =
        0.6931471805599453094172321214581765680755; // +
    const double sqrt_of_ln_of_2 =
        0.8325546111576977563531646448952010476306; // +
    const double pi =
        3.1415926535897932384626433832795028841972; // -
    const double sqrt_of_pi =
        1.7724538509055160272981674833411451827975; // +

    // Functions
    double gauss(double x, double b);
    double lorentz(double x, double d);
    double pseudovoigt(double x, double b, double d, double z);
    double voigt(double x, double b, double d);

    // Function-like classes
    class gaumm_pf;
    class gauss_pf;
    class test;
    class voigt_pf;

    // Funktion-like class templates
    template<class profile_function> class superposition;
}

inline
double
RQ::gauss(double x, double b)
{
    using std::exp;

    return (1.0 / (sqrt_of_pi * b)) * exp(-(x / b) * (x / b));
}

inline
double
RQ::lorentz(double x, double b)
{
    return 1.0 / ((pi * b) * (1.0 + (x / b) * (x / b)));
}

inline
double
RQ::pseudovoigt(double x, double b, double d, double z)
{
    return (1.0 - z) * gauss(x, b) + z * lorentz(x, d);
}

class RQ::gaumm_pf {
public:
    static const size_t parameters = 8;

    gaumm_pf();
    gaumm_pf(const double a[]);
        // a[0] laboratory wavelength (Angstrom)
        // a[1] oscillator strength
        // a[2] cosmological redshift
        // a[3] radial velocity (kilometer per second)
        // a[4] line broadening velocity (kilometer per second)
        // a[5] particle column density (dex per square centimeter)
        // a[6] relativistic correction coefficient
        // a[7] variability of the fine-structure constant (1.0e-05)
    ~gaumm_pf();

    double operator()(double x) const;
    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double
RQ::gaumm_pf::operator()(double x) const
{
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * gauss(x - y, b) : 0.0;
}

inline
double
RQ::gaumm_pf::center() const
{
    return y;
}

class RQ::gauss_pf {
public:
    static const size_t parameters = 6;

    gauss_pf();
    gauss_pf(const double a[]);
        // a[0] laboratory wavelength (Angstrom)
        // a[1] oscillator strength
        // a[2] cosmological redshift
        // a[3] radial velocity (kilometer per second)
        // a[4] line broadening velocity (kilometer per second)
        // a[5] particle column density (dex per square centimeter)
    ~gauss_pf();

    double operator()(double x) const;
    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double
RQ::gauss_pf::operator()(double x) const
{
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * gauss(x - y, b) : 0.0;
}

inline
double
RQ::gauss_pf::center() const
{
    return y;
}

class RQ::test {
public:
    static const size_t parameters = 5;

    test();
    test(const double a[]);
        // a[0] laboratory wavelength (Angstrom)
        // a[1] oscillator strength
        // a[2] observed wavelength (Angstrom)
        // a[3] line broadening velocity (kilometer per second)
        // a[4] particle column density (dex per square centimeter)
    ~test();

    double operator()(double x) const;
    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double
RQ::test::operator()(double x) const
{
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * gauss(x - y, b) : 0.0;
}

inline
double
RQ::test::center() const
{
    return y;
}

class RQ::voigt_pf {
public:
    static const size_t parameters = 7;

    voigt_pf();
    voigt_pf(const double a[]);
        // a[0] laboratory wavelength (Angstrom)
        // a[1] oscillator strength
        // a[2] cosmological redshift
        // a[3] radial velocity (kilometer per second)
        // a[4] line broadening velocity (kilometer per second)
        // a[5] particle column density (dex per square centimeter)
        // a[6] damping constant (Hertz)
    ~voigt_pf();

    double operator()(double x) const;
    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double d; // Lorentzian width (Angstrom)
    double y; // central wavelength (Angstrom)
    double z; // pseudo-voigt function mixing parameter
};

inline
double
RQ::voigt_pf::operator()(double x) const
{
    return c * pseudovoigt(x - y, b, d, z);
}

inline
double
RQ::voigt_pf::center() const
{
    return y;
}

template<class profile_function>
class RQ::superposition {
public:
    superposition();
    superposition(size_t n, const double a[]);
    ~superposition();

    double operator()(double x) const;
    void assign(size_t n, const double a[]);

private:
    std::vector<profile_function> p;
};

template <class profile_function>
RQ::superposition<profile_function>::superposition()
    :   p()
{
}

template <class profile_function>
RQ::superposition<profile_function>::superposition(size_t n, const double a[])
    :   p(n)
{
    for (size_t i = 0; i < n; ++i, a += profile_function::parameters)
        p[i].assign(a);
}

template <class profile_function>
RQ::superposition<profile_function>::~superposition()
{
}

template<class profile_function>
double
RQ::superposition<profile_function>::operator()(double x) const
{
    double d = 0.0;

    for (size_t i = 0; i < p.size(); ++i)
        d += p[i](x);

    return d;
}

template <class profile_function>
void
RQ::superposition<profile_function>::assign(size_t n, const double a[])
{
    this->resize(n);

    for (size_t i = 0; i < n; ++i, a += profile_function::parameters)
        p[i].assign(a);
}

#endif // RQ_PROFILE_FUNCTIONS_H

// References
//
// T. Ida, M. Ando, H. Toraya (2000)
//   Extended pseudo-Voigt function for approximating the Voigt profile,
//   J. Appl. Chryst., 33, 1311, ISSN 0021-8898
