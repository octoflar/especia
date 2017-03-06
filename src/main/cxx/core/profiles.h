/// @file profiles.h
/// Profile funtions.
/// Copyright (c) 2016 Ralf Quast
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
         * @param[in] b The width of the Gaussian (arbitrary unit).
         * @param[in] d The width of the Lorentzian (arbitrary unit).
         */
        Pseudo_Voigt(const Real_t &b = 1.0, const Real_t &d = 1.0);

        /**
         * Destructor.
         */
        ~Pseudo_Voigt();

        /**
         * Returns the value of the pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value (arbitrary unit).
         * @return the value of the pseudo-Voigt approximation at @c x.
         */
        Real_t operator()(const Real_t &x) const;

    private:
        const Real_t u;
        const Real_t r;
        const Real_t gamma_g;
        const Real_t gamma_l;
        const Real_t eta;

        static const Real_t c_g;
        static const Real_t c_l;
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
         * @param[in] b The width of the Gaussian (arbitrary unit).
         * @param[in] d The width of the Lorentzian (arbitrary unit).
         */
        Extended_Pseudo_Voigt(const Real_t &b = 1.0, const Real_t &d = 1.0);

        /**
         * Destructor.
         */
        ~Extended_Pseudo_Voigt();

        /**
         * Returns the value of the extended pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value (arbitrary unit).
         * @return the value of the extended pseudo-Voigt approximation at @c x.
         */
        Real_t operator()(const Real_t &x) const;

    private:
        const Real_t u;
        const Real_t r;
        const Real_t gamma_g;
        const Real_t gamma_l;
        const Real_t gamma_i;
        const Real_t gamma_p;
        const Real_t eta_l;
        const Real_t eta_i;
        const Real_t eta_p;

        static const Real_t c_g;
        static const Real_t c_l;
        static const Real_t c_i;
        static const Real_t c_p;
    };

    /**
     * The (Doppler) profile to infer the variation of the fine-structure constant
     * alpha by means of a many-multiplet analysis.
     *
     * Further reading:
     *
     * R. Quast, D. Reimers and S. A. Levshakov (2004).
     *   *Probing the variability of the fine-structure constant with the VLT/UVES.*
     *   Astronomy and Astrophysics, 415, L7.
     *   http://dx.doi.org/10.1051/0004-6361:20040013
     */
    class Many_Multiplet {
    public:
        /**
         * The number of parameters.
         */
        static const Nnum_t parameter_count = 8;

        /**
         * Default constructor.
         */
        Many_Multiplet();

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
         * @c q[7] The variation of the fine-structure constant (1E-6)
         * @endparblock
         */
        Many_Multiplet(const Real_t q[]);

        /**
         * Destructor.
         */
        ~Many_Multiplet();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return the value of the Doppler profile at @c x.
         */
        Real_t operator()(const Real_t &x) const;

    private:
        /**
         * The modified rest wavelength (Angstrom).
         */
        const Real_t u;

        /**
         * The central wavelength (Angstrom).
         */
        const Real_t c;

        /**
         * The Doppler width (Angstrom).
         */
        const Real_t b;

        /**
         * The amplitude.
         */
        const Real_t a;

        static const Real_t C0;
        static const Real_t C1;
    };


    /**
     * The Doppler profile to model intergalactic absorption lines.
     */
    class Intergalactic_Doppler {
    public:
        /**
         * The number of parameters.
         */
        static const Nnum_t parameter_count = 6;

        /**
         * Default constructor.
         */
        Intergalactic_Doppler();

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
        Intergalactic_Doppler(const Real_t q[]);

        /**
         * Destructor.
         */
        ~Intergalactic_Doppler();

        /**
         * Returns the value of the Doppler profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return the value of the Doppler profile at @c x.
         */
        Real_t operator()(const Real_t &x) const;

    private:
        /**
         * The central wavelength (Angstrom).
         */
        const Real_t c;

        /**
         * The Doppler width (Angstrom).
         */
        const Real_t b;

        /**
         * The amplitude.
         */
        const Real_t a;

        static const Real_t C0;
        static const Real_t C1;
    };


    /**
     * The Voigt profile to model intergalactic spectral lines.
     *
     * @tparam A The strategy to approximate the Voigt function.
     */
    template<class A>
    class Intergalactic_Voigt {
    public:
        /**
         * The number of parameters.
         */
        static const Nnum_t parameter_count = 7;

        /**
         * Default constructor.
         */
        Intergalactic_Voigt()
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
        Intergalactic_Voigt(const Real_t q[])
                : c(q[0] * (1.0 + q[2]) * (1.0 + q[3] / C0)),
                  a(C1 * q[1] * pow(10.0, q[5]) * (q[0] * c)),
                  approximation(q[4] * c / C0, C2 * q[6] * (q[0] * c)) {
        }

        /**
         * Destructor.
         */
        ~Intergalactic_Voigt() {
        }

        /**
         * Returns the value of the Voigt profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom)
         * @return the value of the Voigt profile at @c x.
         */
        Real_t operator()(const Real_t &x) const {
            return a * approximation(x - c);
        };

    private:
        /**
         * The central wavelength (Angstrom).
         */
        const Real_t c;

        /**
         * The amplitude.
         */
        const Real_t a;

        /**
         * The approximation.
         */
        const A approximation;

        static const Real_t C0;
        static const Real_t C1;
        static const Real_t C2;
    };

    template<class A>
    const Real_t Intergalactic_Voigt<A>::C0 = 1.0E-03 * speed_of_light;

    template<class A>
    const Real_t Intergalactic_Voigt<A>::C1 = 1.0E-06 * sqr(elementary_charge) /
                                              (4.0 * electric_constant * electron_mass * sqr(speed_of_light));

    template<class A>
    const Real_t Intergalactic_Voigt<A>::C2 = 1.0E-10 / (4.0 * pi * speed_of_light);


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
        Superposition(Nnum_t n, const Real_t q[])
                : profiles() {
            profiles.reserve(n);
            for (Nnum_t i = 0; i < n; ++i, q += P::parameter_count) {
                profiles.push_back(P(q));
            }
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
         * @return the value of the profile superposition at @c x.
         */
        Real_t operator()(const Real_t &x) const {
            Real_t d = 0.0;

            for (Nnum_t i = 0; i < profiles.size(); ++i) {
                d += profiles[i](x);
            }

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
     * @param[in] f The function.
     * @param[in] x The wavelength relative to the center of the profile.
     * @param[in] b The width of the profile.
     * @param[in] c The truncation parameter.
     * @return the value of the profile function at @c x, if the absolute value
     *         of @c x is less than @c c widths, zero otherwise.
     */
    template<class F>
    inline
    Real_t truncate(const F &f, const Real_t &x, const Real_t &b, const Real_t &c) {
        using std::abs;

        return abs(x) < c * b ? f(x, b) : 0.0;
    }
}

#endif // ESPECIA_PROFILES_H
