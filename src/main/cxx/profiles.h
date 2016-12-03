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
#ifndef ESPECIA_PROFILES_H
#define ESPECIA_PROFILES_H

#include <cmath>
#include <cstddef>
#include <vector>
#include "base.h"

namespace especia {

    // Function-like classes
    class doppler_ig; // Doppler profile
    class doppler_mm; // Doppler profile for many-multiplets analysis
    class doppler_is; // Doppler profile for interstellar lines
    class voigt_ig;   // Voigt profile

    // Function-like class templates

    template<class profile_function>
    class superposition;
}


class especia::doppler_mm {
public:
    static const size_t parameters = 8;

    doppler_mm();

    doppler_mm(const double a[]);

    // a[0] laboratory wavelength (Angstrom)
    // a[1] oscillator strength
    // a[2] cosmological redshift
    // a[3] radial velocity (km s-1)
    // a[4] line broadening velocity (km s-1)
    // a[5] decadic logarithm of the particle column number density (cm-2)
    // a[6] relativistic correction coefficient
    // a[7] variability of the fine-structure constant (1.0e-05)
    ~doppler_mm();

    double operator()(double x) const;

    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double especia::doppler_mm::center() const {
    return y;
}

class especia::doppler_ig {
public:
    static const size_t parameters = 6;

    doppler_ig();

    doppler_ig(const double a[]);

    // a[0] laboratory wavelength (Angstrom)
    // a[1] oscillator strength
    // a[2] cosmological redshift
    // a[3] radial velocity (km s-1)
    // a[4] line broadening velocity (km s-1)
    // a[5] decadic logarithm of the particle column number density (cm-2)
    ~doppler_ig();

    double operator()(double x) const;

    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double especia::doppler_ig::center() const {
    return y;
}

class especia::doppler_is {
public:
    static const size_t parameters = 5;

    doppler_is();

    doppler_is(const double a[]);

    // a[0] laboratory wavelength (Angstrom)
    // a[1] oscillator strength
    // a[2] radial velocity (km s-1)
    // a[3] line broadening velocity (km s-1)
    // a[4] decadic logarithm of the particle column number density (cm-2)
    ~doppler_is();

    double operator()(double x) const;

    double center() const;

    void assign(const double a[]);

private:
    double b; // Doppler width (Angstrom)
    double c; // amplitude
    double y; // central wavelength (Angstrom)
};

inline
double especia::doppler_is::center() const {
    return y;
}

class especia::voigt_ig {
public:
    static const size_t parameters = 7;

    voigt_ig();

    voigt_ig(const double a[]);

    // a[0] laboratory wavelength (Angstrom)
    // a[1] oscillator strength
    // a[2] cosmological redshift
    // a[3] radial velocity (km s-1)
    // a[4] line broadening velocity (km s-1)
    // a[5] decadic logarithm of the particle column number density (cm-2)
    // a[6] damping constant (s-1)
    ~voigt_ig();

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
double especia::voigt_ig::center() const {
    return y;
}

template<class profile_function>
class especia::superposition {
public:
    superposition()
            : p() {
    }

    superposition(size_t n, const double a[])
            : p(n) {
        for (size_t i = 0; i < n; ++i, a += profile_function::parameters)
            p[i].assign(a);
    }

    ~superposition() {
    }

    double operator()(double x) const {
        double d = 0.0;

        for (size_t i = 0; i < p.size(); ++i)
            d += p[i](x);

        return d;
    }

    void assign(size_t n, const double a[]) {
        p.resize(n);

        for (size_t i = 0; i < n; ++i, a += profile_function::parameters)
            p[i].assign(a);
    }

private:
    std::vector<profile_function> p;
};

#endif // ESPECIA_PROFILES_H

// References
//
// T. Ida, M. Ando, H. Toraya (2000)
//   Extended pseudo-Voigt function for approximating the Voigt profile,
//   J. Appl. Chryst., 33, 1311, ISSN 0021-8898
