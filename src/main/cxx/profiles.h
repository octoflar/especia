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

    class pseudo_voigt;
    class extended_pseudo_voigt;

    // Function-like classes
    class doppler_ig; // Doppler profile
    class doppler_mm; // Doppler profile for many-multiplets analysis
    class doppler_is; // Doppler profile for interstellar lines

    template <class approximation>
    class voigt_ig;   // Voigt profile

    // Function-like class templates

    template<class profile_function>
    class superposition;
}


class especia::pseudo_voigt {
public:
    pseudo_voigt(double b = 1.0, double d = 1.0);

    ~pseudo_voigt();

    double operator()(const double &x) const;

private:
    double gamma_g;
    double gamma_l;
    double eta;

    static const double c_g;
    static const double c_l;
};


class especia::extended_pseudo_voigt {
public:
    extended_pseudo_voigt(double b = 1.0, double d = 1.0);

    ~extended_pseudo_voigt();

    double operator()(const double &x) const;

private:
    double gamma_g;
    double gamma_l;
    double gamma_i;
    double gamma_p;
    double eta_l;
    double eta_i;
    double eta_p;

    static const double c_g;
    static const double c_l;
    static const double c_i;
    static const double c_p;
};


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

    double operator()(double x) const;;

    void assign(const double a[]);

private:
    double y; // central wavelength (Angstrom)
    double b; // Doppler width (Angstrom)
    double c; // amplitude
};


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

    double operator()(double x) const;;

    void assign(const double a[]);

private:
    double y; // central wavelength (Angstrom)
    double b; // Doppler width (Angstrom)
    double c; // amplitude
};


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

    double operator()(double x) const;;

    void assign(const double a[]);

private:
    double y; // central wavelength (Angstrom)
    double b; // Doppler width (Angstrom)
    double c; // amplitude
};


template<class approximation>
class especia::voigt_ig {
public:
    static const size_t parameters = 7;

    voigt_ig() : c(0.0), y(0.0), f(1.0, 1.0) {
    };

    // a[0] laboratory wavelength (Angstrom)
    // a[1] oscillator strength
    // a[2] cosmological redshift
    // a[3] radial velocity (km s-1)
    // a[4] line broadening velocity (km s-1)
    // a[5] decadic logarithm of the particle column number density (cm-2)
    // a[6] damping constant (s-1)
    voigt_ig(const double a[]){
        assign(a);
    }

    ~voigt_ig(){
    }

    double operator()(double x) const {
        return c * f(x - y);
    };

    void assign(const double a[]) {
        y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / SPEED_OF_LIGHT);
        c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);

        const double b = a[4] * y / SPEED_OF_LIGHT;
        const double d = 2.65442e-20 * a[6] * (a[0] * y);

        f = approximation(b, d);
    }

private:
    double y; // central wavelength (Angstrom)
    double c; // amplitude
    approximation f;
};


template<class profile>
class especia::superposition {
public:
    superposition()
            : p() {
    }

    superposition(size_t n, const double a[])
            : p(n) {
        for (size_t i = 0; i < n; ++i, a += profile::parameters)
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

        for (size_t i = 0; i < n; ++i, a += profile::parameters)
            p[i].assign(a);
    }

private:
    std::vector<profile> p;
};

#endif // ESPECIA_PROFILES_H

// References
//
// T. Ida, M. Ando, H. Toraya (2000)
//   Extended pseudo-Voigt function for approximating the Voigt profile,
//   J. Appl. Chryst., 33, 1311, ISSN 0021-8898
