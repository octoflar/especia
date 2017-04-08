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
        Pseudo_Voigt(const R_type &b = 0.5, const R_type &d = 0.5);

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
        R_type operator()(const R_type &x) const;

    private:
        const R_type u;
        const R_type r;
        const R_type gamma_g;
        const R_type gamma_l;
        const R_type eta;

        static const R_type c_g;
        static const R_type c_l;
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
        Extended_Pseudo_Voigt(const R_type &b = 0.5, const R_type &d = 0.5);

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
        R_type operator()(const R_type &x) const;

    private:
        const R_type u;
        const R_type r;
        const R_type gamma_g;
        const R_type gamma_l;
        const R_type gamma_i;
        const R_type gamma_p;
        const R_type eta_l;
        const R_type eta_i;
        const R_type eta_p;

        static const R_type c_g;
        static const R_type c_l;
        static const R_type c_i;
        static const R_type c_p;
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
        Many_Multiplet(const R_type q[]);

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
        R_type operator()(const R_type &x) const;

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        R_type get_center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        R_type get_redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static N_type get_parameter_count() {
            return parameter_count;
        };

    private:
        /**
         * The modified rest wavelength (Angstrom).
         */
        const R_type u;

        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const R_type z;

        /**
         * The central wavelength (Angstrom).
         */
        const R_type c;

        /**
         * The Doppler width (Angstrom).
         */
        const R_type b;

        /**
         * The amplitude.
         */
        const R_type a;

        /**
         * The number of parameters.
         */
        static const N_type parameter_count = 8;

        static const R_type c0;
        static const R_type c1;
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
        Intergalactic_Doppler(const R_type q[]);

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
        R_type operator()(const R_type &x) const;

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        R_type get_center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        R_type get_redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static N_type get_parameter_count() {
            return parameter_count;
        };

    private:
        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const R_type z;

        /**
         * The central wavelength (Angstrom).
         */
        const R_type c;

        /**
         * The Doppler width (Angstrom).
         */
        const R_type b;

        /**
         * The amplitude.
         */
        const R_type a;

        /**
         * The number of parameters.
         */
        static const N_type parameter_count = 6;

        static const R_type c0;
        static const R_type c1;
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
                : z(1.0), a(1.0), c(0.0), approximation(1.0, 1.0) {
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
        Intergalactic_Voigt(const R_type q[])
                : z((1.0 + q[2]) * (1.0 + q[3] / c0)),
                  c(q[0] * z),
                  a(c1 * q[1] * pow(10.0, q[5]) * (q[0] * c)),
                  approximation(q[4] * c / c0, c2 * q[6] * (q[0] * c)) {
        }

        /**
         * The destructor.
         */
        ~Intergalactic_Voigt() {
        }

        /**
         * Returns the optical depth of the profile at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile at @c x.
         */
        R_type operator()(const R_type &x) const {
            return a * approximation(x - c);
        };

        /**
         * Returns the central wavelength of the profile.
         *
         * @return the central wavelength of the profile (Angstrom).
         */
        R_type get_center() const {
            return c;
        }

        /**
         * Returns the redshift factor of the profile due to cosmology and proper motion.
         *
         * @return the redshift factor.
         */
        R_type get_redshift_factor() const {
            return z;
        }

        /**
         * Returns the number of parameters.
         */
        static N_type get_parameter_count() {
            return parameter_count;
        };

    private:
        /**
         * The redshift factor due to cosmology and proper motion.
         */
        const R_type z;

        /**
         * The central wavelength (Angstrom).
         */
        const R_type c;

        /**
         * The amplitude.
         */
        const R_type a;

        /**
         * The approximation.
         */
        const A approximation;

        /**
         * The number of parameters.
         */
        static const N_type parameter_count = 7;

        static const R_type c0;
        static const R_type c1;
        static const R_type c2;
    };

    template<class A>
    const R_type Intergalactic_Voigt<A>::c0 = 1.0E-03 * speed_of_light;

    template<class A>
    const R_type Intergalactic_Voigt<A>::c1 = 1.0E-06 * sq(elementary_charge) /
                                              (4.0 * electric_constant * electron_mass * sq(speed_of_light));

    template<class A>
    const R_type Intergalactic_Voigt<A>::c2 = 1.0E-10 / (4.0 * pi * speed_of_light);


    /**
     * The superposition of many optical depth profiles.
     *
     * @tparam T The profile type.
     *
     * @remark This class is thread safe, if the profile type is thead safe.
     */
    template<class T>
    class Superposition {
    public:
        /**
         * Constructs a new superposition of profiles with the parameter values specified.
         *
         * @param[in] n The number of profiles.
         * @param[in] q The vector of parameter values. The semantics of parameter values and the
         *              number of parameters per component are defined by the profile type.
         */
        Superposition(N_type n, const R_type q[])
                : profiles() {
            profiles.reserve(n);
            for (N_type i = 0; i < n; ++i, q += T::get_parameter_count()) {
                profiles.push_back(T(q));
            }
        }

        /**
         * The destructor.
         */
        ~Superposition() {
        }

        /**
         * Returns the optical depth of the profile superposition at a given wavelength.
         *
         * @param[in] x The wavelength (Angstrom).
         * @return the optical depth of the profile superposition at @c x.
         */
        R_type operator()(const R_type &x) const {
            R_type d = 0.0;

            for (N_type i = 0; i < profiles.size(); ++i) {
                d += profiles[i](x);
            }

            return d;
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
        Equivalent_Width_Calculator(const Integrate &integrator) : integrator(integrator) {

        }

        /**
         * The destructor.
         */
        ~Equivalent_Width_Calculator() {

        }

        /**
         * Calculates the equivalent width of an optical depth profile.
         *
         * @tparam T The profile type.
         *
         * @param t The profile.
         * @return the equivalent width (milli Angstrom).
         */
        template<class T>
        R_type calculate(const T &t) const {
            using std::exp;

            const R_type integral = integrator.integrate_semi_infinite(
                    [&t](R_type x) -> R_type { return kilo * (1.0 - exp(-t(x + t.get_center()))); });

            return 2.0 * integral / t.get_redshift_factor();
        }

    private:
        /**
         * The strategy to integrate the line profile.
         */
        const Integrate integrator;
    };
}

#endif // ESPECIA_PROFILES_H
