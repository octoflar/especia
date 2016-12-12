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
    class Pseudo_Voigt {
    public:
        /**
         * Creates a new pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian.
         * @param[in] d The width of the Lorentzian.
         */
        Pseudo_Voigt(const double &b = 1.0, const double &d = 1.0);

        /**
         * Destructor.
         */
        ~Pseudo_Voigt();

        /**
         * Returns the value of the pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value.
         * @return the value of the pseudo-Voigt approximation at @c x.
         */
        double operator()(const double &x) const;

    private:
        const double u;
        const double r;
        const double gamma_g;
        const double gamma_l;
        const double eta;

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
    class Extended_Pseudo_Voigt {
    public:
        /**
         * Creates a new extended pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian.
         * @param[in] d The width of the Lorentzian.
         */
        Extended_Pseudo_Voigt(const double &b = 1.0, const double &d = 1.0);

        /**
         * Destructor.
         */
        ~Extended_Pseudo_Voigt();

        /**
         * Returns the value of the extended pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value.
         * @return the value of the extended pseudo-Voigt approximation at @c x.
         */
        double operator()(const double &x) const;

    private:
        const double u;
        const double r;
        const double gamma_g;
        const double gamma_l;
        const double gamma_i;
        const double gamma_p;
        const double eta_l;
        const double eta_i;
        const double eta_p;

        static const double c_g;
        static const double c_l;
        static const double c_i;
        static const double c_p;
    };

    /**
     * The Doppler profile to infer the variation of the fine-structure constant
     * alpha by means of a many-multiplet analysis.
     *
     * Further reading:
     *
     * R. Quast, D. Reimers and S. A. Levshakov (2004).
     *   *Probing the variability of the fine-structure constant with the VLT/UVES.*
     *   Astronomy and Astrophysics, 415, L7.
     *   http://dx.doi.org/10.1051/0004-6361:20040013
     */
    class A_Doppler {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 8;

        /**
         * Default constructor.
         */
        A_Doppler();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c q[0] The rest wavelength (Angstrom)
         *
         * @c q[1] The oscillator strength
         *
         * @c q[2] The cosmological redshift
         *
         * @c q[3] The radial velocity (km s-1)
         *
         * @c q[4] The line broadening velocity (km s-1)
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c q[6] The relativistic correction coefficient
         *
         * @c q[7] The variation of the fine-structure constant (1.0e-05)
         * @endparblock
         */
        A_Doppler(const double q[]);

        /**
         * Destructor.
         */
        ~A_Doppler();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Doppler profile at @c x.
         */
        double operator()(double x) const;

    private:
        /**
         * The modified rest wavelength (Angstrom).
         */
        const double u;

        /**
         * The central wavelength (Angstrom).
         */
        const double c;

        /**
         * The Doppler width (Angstrom).
         */
        const double b;

        /**
         * The amplitude.
         */
        const double a;

        static const double C0;
        static const double C1;
    };


    /**
     * The Doppler profile to model intergalactic absorption lines.
     */
    class G_Doppler {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 6;

        /**
         * Default constructor.
         */
        G_Doppler();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c q[0] The rest wavelength (Angstrom)
         *
         * @c q[1] The oscillator strength
         *
         * @c q[2] The cosmological redshift
         *
         * @c q[3] The radial velocity (km s-1)
         *
         * @c q[4] The line broadening velocity (km s-1)
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2)
         * @endparblock
         */
        G_Doppler(const double q[]);

        /**
         * Destructor.
         */
        ~G_Doppler();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Doppler profile at @c x.
         */
        double operator()(double x) const;

    private:
        /**
         * The central wavelength (Angstrom).
         */
        const double c;

        /**
         * The Doppler width (Angstrom).
         */
        const double b;

        /**
         * The amplitude.
         */
        const double a;

        static const double C0;
        static const double C1;
    };


    /**
     * The Voigt profile to model intergalactic spectral lines.
     *
     * @tparam A The strategy to approximate the Voigt function.
     */
    template<class A>
    class G_Voigt {
    public:
        /**
         * The number of parameters.
         */
        static const size_t parameter_count = 7;

        /**
         * Default constructor.
         */
        G_Voigt()
                : a(0.0), c(0.0), approximation(1.0, 1.0) {
        };

        /**
         * Creates a new Voigt profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The vector of parameter values. Its components are:
         *
         * @c q[0] The rest wavelength (Angstrom)
         *
         * @c q[1] The oscillator strength
         *
         * @c q[2] The cosmological redshift
         *
         * @c q[3] The radial velocity (km s-1)
         *
         * @c q[4] The line broadening velocity (km s-1)
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2)
         *
         * @c q[6] The damping constant (s-1)
         * @endparblock
         */
        G_Voigt(const double q[])
                : c(q[0] * (1.0 + q[2]) * (1.0 + q[3] / C0)),
                  a(C1 * q[1] * pow(10.0, q[5]) * (q[0] * c)),
                  approximation(q[4] * c / C0, C2 * q[6] * (q[0] * c)) {
        }

        /**
         * Destructor.
         */
        ~G_Voigt() {
        }

        /**
         * Returns the value of the Voigt profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the Voigt profile at @c x.
         */
        double operator()(double x) const {
            return a * approximation(x - c);
        };

    private:
        /**
         * The central wavelength (Angstrom).
         */
        const double c;

        /**
         * The amplitude.
         */
        const double a;

        /**
         * The approximation.
         */
        const A approximation;

        static const double C0;
        static const double C1;
        static const double C2;
    };

    template<class A>
    const double G_Voigt<A>::C0 = 1.0E-03 * speed_of_light_in_vacuum;

    template<class A>
    const double G_Voigt<A>::C1 = 1.0E-06 * sqr(elementary_charge) /
                                  (4.0 * electric_constant * electron_mass * sqr(speed_of_light_in_vacuum));
    template<class A>
    const double G_Voigt<A>::C2 = 1.0E-10 / (4.0 * pi * speed_of_light_in_vacuum);


    /**
     * The superposition of many profiles.
     *
     * @tparam P The profile type.
     */
    template<class P>
    class Superposition {
    public:
        /**
         * Constructs a new superposition of profiles with the parameter values specified.
         *
         * @param[in] n The number of profiles.
         * @param[in] q The vector of parameter values. The semantics of parameter values and the
         *              number of parameters per component are defined by the profile type.
         */
        Superposition(size_t n, const double q[])
                : profiles() {
            profiles.reserve(n);
            for (size_t i = 0; i < n; ++i, q += P::parameter_count)
                profiles.push_back(P(q));
        }

        /**
         * Destructor.
         */
        ~Superposition() {
        }

        /**
         * Returns the value of the profile superpositon at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return The value of the profile superposition at @c x.
         */
        double operator()(double x) const {
            double d = 0.0;

            for (size_t i = 0; i < profiles.size(); ++i)
                d += profiles[i](x);

            return d;
        }

    private:
        std::vector<P> profiles;
    };

    /**
     * Truncates the support of a given profile function.
     *
     * @tparam F The function type.
     *
     * @param f The function.
     * @param x The wavelength relative to the center of the profile.
     * @param b The width of the profile.
     * @param c The truncation parameter.
     * @return The value of the profile function at @c x, if the
     *         absolute value of @c x is less than @c c, zero
     *         otherwise.
     */
    template<class F>
    inline
    double truncate(const F &f, const double &x, const double& b, const double &c) {
        using std::abs;

        return abs(x) < c * b ? f(x, b) : 0.0;
    }
}

#endif // ESPECIA_PROFILES_H
