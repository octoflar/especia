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

   /**
    * The pseudo-Voigt approximation to the Voigt function. The Voigt function is
    * defined as the convolution of a Gaussian and a Lorentzian function.
    *
    * Further reading:
    *
    * T. Ida, M. Ando, H. Toraya (2000).
    *   *Extended pseudo-Voigt function for approximating the Voigt profile.*
    *   J. Appl. Chryst., 33, 1311, ISSN 0021-8898.
    */
    class pseudo_voigt {
    public:
        /**
         * Creates a new pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian.
         * @param[in] d The width of the Lorentzian.
         */
        pseudo_voigt(double b = 1.0, double d = 1.0);

        /**
         * Destructor.
         */
        ~pseudo_voigt();

        /**
         * Returns the value of the pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value.
         * @return the value of the pseudo-Voigt approximation at @c x.
         */
        double operator()(const double &x) const;

    private:
        double gamma_g;
        double gamma_l;
        double eta;

        static const double c_g;
        static const double c_l;
    };


   /**
    * The extended pseudo-Voigt approximation to the Voigt function. The Voigt function
    * is defined as the convolution of a Gaussian and a Lorentzian function.
    *
    * Further reading:
    *
    * T. Ida, M. Ando, H. Toraya (2000)
    *   *Extended pseudo-Voigt function for approximating the Voigt profile.*
    *   J. Appl. Chryst., 33, 1311, ISSN 0021-8898.
    */
    class extended_pseudo_voigt {
    public:
        /**
         * Creates a new extended pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian.
         * @param[in] d The width of the Lorentzian.
         */
        extended_pseudo_voigt(double b = 1.0, double d = 1.0);

        /**
         * Destructor.
         */
        ~extended_pseudo_voigt();

        /**
         * Returns the value of the extended pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value.
         * @return the value of the extended pseudo-Voigt approximation at @c x.
         */
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

   /**
    * The Doppler profile to infer the variation of the fine-structure constant
    * by means of a many-multiplet analysis.
    *
    * Further reading:
    *
    * R. Quast, D. Reimers and S. A. Levshakov (2004).
    *   *Probing the variability of the fine-structure constant with the VLT/UVES.*
    *   Astronomy and Astrophysics, 415, L7.
    *   doi: http://dx.doi.org/10.1051/0004-6361:20040013
    */
    class doppler_mm {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 8;

        /**
         * Default constructor.
         */
        doppler_mm();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c a[6] The relativistic correction coefficient
         *
         * @c a[7] The variation of the fine-structure constant (1.0e-05)
         * @endparblock
         */
        doppler_mm(const double a[]);

        /**
         * Destructor.
         */
        ~doppler_mm();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Doppler profile at @c x.
         */
        double operator()(double x) const;

        /**
         * Assigns a new set of parameter values to this Doppler profile.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c a[6] The relativistic correction coefficient
         *
         * @c a[7] The variation of the fine-structure constant (1.0e-05)
         * @endparblock
         */
        void assign(const double a[]);

    private:
        double y; // central wavelength (Angstrom)
        double b; // Doppler width (Angstrom)
        double c; // amplitude
    };


   /**
    * The Doppler profile to model intergalactic absorption lines.
    */
    class doppler_ig {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 6;

        /**
         * Default constructor.
         */
        doppler_ig();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         * @endparblock
         */
        doppler_ig(const double a[]);

        /**
         * Destructor.
         */
        ~doppler_ig();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Doppler profile at @c x.
         */
        double operator()(double x) const;

        /**
         * Assigns a new set of parameter values to this Doppler profile.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         * @endparblock
         */
        void assign(const double a[]);

    private:
        double y; // central wavelength (Angstrom)
        double b; // Doppler width (Angstrom)
        double c; // amplitude
    };


   /**
    * The Voigt profile to model intergalactic spectral lines.
    *
    * @tparam approximation The approximation to the Voigt function.
    */
    template<class approximation>
    class voigt_ig {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 7;

        /**
         * Default constructor.
         */
        voigt_ig() : c(0.0), y(0.0), f(1.0, 1.0) {
        };

        /**
         * Creates a new Voigt profile with the parameter values specified.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c a[6] The damping constant (s-1)
         * @endparblock
         */
        voigt_ig(const double a[]) {
            assign(a);
        }

        /**
         * Destructor.
         */
        ~voigt_ig() {
        }

        /**
         * Returns the value of the Voigt profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Voigt profile at @c x.
         */
        double operator()(double x) const {
            return c * f(x - y);
        };

        /**
         * Assigns a new set of parameter values to this Voigt profile.
         *
         * @param[in] a
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c a[0] The laboratory wavelength (Angstrom)
         *
         * @c a[1] The oscillator strength
         *
         * @c a[2] The cosmological redshift
         *
         * @c a[3] The radial velocity (km s-1)
         *
         * @c a[4] The line broadening velocity (km s-1)
         *
         * @c a[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c a[6] The damping constant (s-1)
         * @endparblock
         */
        void assign(const double a[]) {
            y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
            c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);

            const double b = a[4] * y / speed_of_light;
            const double d = 2.65442e-20 * a[6] * (a[0] * y);

            f = approximation(b, d);
        }

    private:
        double y; // central wavelength (Angstrom)
        double c; // amplitude
        approximation f;
    };


   /**
    * The superposition of many profiles.
    *
    * @tparam profile The profile type.
    */
    template<class profile>
    class superposition {
    public:
        /**
         * Default constructor.
         */
        superposition()
                : p() {
        }

        /**
         * Constructs a new superposition of profiles with the parameter values specified.
         *
         * @param[in] n The number of profiles.
         * @param[in] a The vector of parameter values. The semantics of parameter values and the
         *              number of parameters per component are defined by the profile type.
         */
        superposition(size_t n, const double a[])
                : p(n) {
            for (size_t i = 0; i < n; ++i, a += profile::parameter_count)
                p[i].assign(a);
        }

        /**
         * Destructor.
         */
        ~superposition() {
        }

        /**
         * Returns the value of the profile superpositon at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the profile superposition at @c x.
         */
        double operator()(double x) const {
            double d = 0.0;

            for (size_t i = 0; i < p.size(); ++i)
                d += p[i](x);

            return d;
        }

        /**
         * Assigns a new set of parameter values to this superposition.
         *
         * @param[in] n The number of profiles.
         * @param[in] a The vector of parameter values. The semantics of parameter values and the
         *              number of parameters per component are defined by the profile type.
         */
        void assign(size_t n, const double a[]) {
            p.resize(n);

            for (size_t i = 0; i < n; ++i, a += profile::parameter_count)
                p[i].assign(a);
        }

    private:
        std::vector<profile> p;
    };

}

#endif // ESPECIA_PROFILES_H
