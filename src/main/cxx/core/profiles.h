/// @file profiles.h
/// Profile funtions.
/// Copyright (c) 2017 Ralf Quast
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
     *
     * @remark This class is thread safe.
     */
    class Pseudo_Voigt {
    public:
        /**
         * Creates a new pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian (arbitrary unit).
         * @param[in] d The width of the Lorentzian (arbitrary unit).
         */
        explicit Pseudo_Voigt(const real &b = 0.5, const real &d = 0.5);

        /**
         * The destructor.
         */
        ~Pseudo_Voigt();

        /**
         * Returns the value of the pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value (arbitrary unit).
         * @return the value of the pseudo-Voigt approximation at @c x.
         */
        real operator()(const real &x) const;

    private:
        const real u;
        const real r;
        const real gamma_g;
        const real gamma_l;
        const real eta;

        static const real c_g;
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
     *
     * @remark This class is thread safe.
     */
    class Extended_Pseudo_Voigt {
    public:
        /**
         * Creates a new extended pseudo-Voigt approximation to the Voigt function.
         *
         * @param[in] b The width of the Gaussian (arbitrary unit).
         * @param[in] d The width of the Lorentzian (arbitrary unit).
         */
        explicit Extended_Pseudo_Voigt(const real &b = 0.5, const real &d = 0.5);

        /**
         * The destructor.
         */
        ~Extended_Pseudo_Voigt();

        /**
         * Returns the value of the extended pseudo-Voigt approximation at a given abscissa value.
         *
         * @param[in] x The abscissa value (arbitrary unit).
         * @return the value of the extended pseudo-Voigt approximation at @c x.
         */
        real operator()(const real &x) const;

    private:
        const real u;
        const real r;
        const real gamma_g;
        const real gamma_l;
        const real gamma_i;
        const real gamma_p;
        const real eta_l;
        const real eta_i;
        const real eta_p;

        static const real c_g;
        static const real c_i;
        static const real c_p;
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
     *   https://dx.doi.org/10.1051/0004-6361:20040013
     *
     * @remark This class is thread safe.
     */
    class Many_Multiplet {
    public:
        /**
         * Default constructor.
         */
        Many_Multiplet();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The parameter values. The parameters are:
         *
         * @c q[0] The rest wavelength (Angstrom).
         *
         * @c q[1] The oscillator strength.
         *
         * @c q[2] The cosmological redshift.
         *
         * @c q[3] The radial velocity (km s-1).
         *
         * @c q[4] The line broadening velocity (km s-1).
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2).
         *
         * @c q[6] The relativistic correction coefficient.
         *
         * @c q[7] The variation of the fine-structure constant (1E-6).
         * @endparblock
         */
        explicit Many_Multiplet(const real q[]);

        /**
         * The destructor.
         */
        ~Many_Multiplet();

        /**
         * Returns the optical depth of the profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile at @c x.
         */
        real operator()(const real &x) const;

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        real center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        real redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static natural parameter_count() {
            return n;
        };

    private:
        /**
         * The modified rest wavelength (Angstrom).
         */
        const real u;

        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const real z;

        /**
         * The central wavelength (Angstrom).
         */
        const real c;

        /**
         * The Doppler width (Angstrom).
         */
        const real b;

        /**
         * The amplitude.
         */
        const real a;

        /**
         * The number of parameters.
         */
        static const natural n = 8;

        static const real c0;
        static const real c1;
    };


    /**
     * The Doppler profile to model intergalactic absorption lines.
     *
     * @remark This class is thread safe.
     */
    class Intergalactic_Doppler {
    public:
        /**
         * Default constructor.
         */
        Intergalactic_Doppler();

        /**
         * Creates a new Doppler profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The parameter values. The parameters are:
         *
         * @c q[0] The rest wavelength (Angstrom).
         *
         * @c q[1] The oscillator strength.
         *
         * @c q[2] The cosmological redshift.
         *
         * @c q[3] The radial velocity (km s-1).
         *
         * @c q[4] The line broadening velocity (km s-1).
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2).
         * @endparblock
         */
        explicit Intergalactic_Doppler(const real q[]);

        /**
         * The destructor.
         */
        ~Intergalactic_Doppler();

        /**
         * Returns the optical depth of the profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile at @c x.
         */
        real operator()(const real &x) const;

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        real center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        real redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static natural parameter_count() {
            return n;
        };

    private:
        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const real z;

        /**
         * The central wavelength (Angstrom).
         */
        const real c;

        /**
         * The Doppler width (Angstrom).
         */
        const real b;

        /**
         * The amplitude.
         */
        const real a;

        /**
         * The number of parameters.
         */
        static const natural n = 6;

        static const real c0;
        static const real c1;
    };


    /**
     * The Voigt profile to model intergalactic spectral lines.
     *
     * @tparam A The strategy to approximate the Voigt function.
     *
     * @remark This class is thread safe.
     */
    template<class A>
    class Intergalactic_Voigt {
    public:
        /**
         * Default constructor.
         */
        Intergalactic_Voigt()
                : z(1.0), c(0.0), a(1.0), approximation(1.0, 1.0) {
        };

        /**
         * Creates a new Voigt profile with the parameter values specified.
         *
         * @param[in] q
         * @parblock
         * The parameter values. The parameters are:
         *
         * @c q[0] The rest wavelength (Angstrom).
         *
         * @c q[1] The oscillator strength.
         *
         * @c q[2] The cosmological redshift.
         *
         * @c q[3] The radial velocity (km s-1).
         *
         * @c q[4] The line broadening velocity (km s-1).
         *
         * @c q[5] The decadic logarithm of the particle column number density (cm-2).
         *
         * @c q[6] The damping constant (s-1).
         * @endparblock
         */
        explicit Intergalactic_Voigt(const real q[])
                : z((1.0 + q[2]) * (1.0 + q[3] / c0)),
                  c(q[0] * z),
                  a(c1 * q[1] * std::pow(10.0, q[5]) * (q[0] * c)),
                  approximation(q[4] * c / c0, c2 * q[6] * (q[0] * c)) {
        }

        /**
         * The destructor.
         */
        ~Intergalactic_Voigt() = default;

        /**
         * Returns the optical depth of the profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile at @c x.
         */
        real operator()(const real &x) const {
            return a * approximation(x - c);
        };

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        real center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        real redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static natural parameter_count() {
            return n;
        };

    private:
        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const real z;

        /**
         * The central wavelength (Angstrom).
         */
        const real c;

        /**
         * The amplitude.
         */
        const real a;

        /**
         * The approximation.
         */
        const A approximation;

        /**
         * The number of parameters.
         */
        static const natural n = 7;

        static const real c0;
        static const real c1;
        static const real c2;
    };

    template<class A>
    const real Intergalactic_Voigt<A>::c0 = 1.0E-03 * speed_of_light;

    template<class A>
    const real Intergalactic_Voigt<A>::c1 = 1.0E-06 * sq(elementary_charge) / // NOLINT
                                                      (4.0 * electric_constant * electron_mass * sq(speed_of_light));

    template<class A>
    const real Intergalactic_Voigt<A>::c2 = 1.0E-10 / (4.0 * pi * speed_of_light);


    /**
     * The superposition of many optical depth profiles.
     *
     * @tparam T The profile type.
     *
     * @remark This class is thread safe, if the profile type is thread safe.
     */
    template<class T>
    class Superposition {
    public:
        /**
         * Constructs a new superposition of profiles with the parameter values specified.
         *
         * @param[in] n The number of profiles.
         * @param[in] q The parameter values. The semantics of parameter values and the
         * number of parameters per component are defined by the profile type.
         */
        Superposition(natural n, const real q[]) : profiles() {
            profiles.reserve(n);
            for (natural i = 0; i < n; ++i, q += T::parameter_count()) {
                profiles.emplace_back(q);
            }
        }

        /**
         * The destructor.
         */
        ~Superposition() = default;

        /**
         * Returns the optical depth of the profile superposition at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile superposition at @c x.
         */
        real operator()(const real &x) const {
            real t = 0.0;

            for (const T &profile : profiles) {
                t += profile(x);
            }

            return t;
        }

    private:
        /**
         * The line profiles.
         */
        std::vector<T> profiles;
    };


    /**
     * Calculates the equivalent width of an optical depth profile.
     *
     * @tparam Integrate The strategy to integrate the line profile.
     */
    template<class Integrate>
    class Equivalent_Width_Calculator {
    public:
        /**
         * The default constructor.
         */
        Equivalent_Width_Calculator() : integrator(Integrate()) {

        }

        /**
         * Constructs a new instance of this class using the integrator supplied as argument.
         *
         * @param integrator The integrator.
         */
        explicit Equivalent_Width_Calculator(const Integrate &integrator) : integrator(integrator) {

        }

        /**
         * The destructor.
         */
        ~Equivalent_Width_Calculator() = default;

        /**
         * Calculates the rest equivalent width of an optical depth profile.
         *
         * @tparam Function The profile function type.
         *
         * @param f The profile function.
         * @return the equivalent width (milli Angstrom).
         */
        template<class Function>
        real calculate(const Function &f) const {
            using std::exp;

            const real integral = integrator.integrate_semi_infinite(
                    [&f](real x) -> real { return kilo * (1.0 - exp(-f(x + f.center()))); });

            return 2.0 * integral / f.redshift_factor();
        }

    private:
        /**
         * The strategy to integrate the line profile.
         */
        const Integrate integrator;
    };
}

#endif // ESPECIA_PROFILES_H

