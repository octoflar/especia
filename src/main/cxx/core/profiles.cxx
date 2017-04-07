/// @file profiles.cxx
/// Profile functions.
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
#include "profiles.h"

using std::abs;
using std::cosh;
using std::exp;
using std::log;
using std::pow;
using std::sqrt;

using especia::N_type;
using especia::R_type;
using especia::micro;
using especia::pi;
using especia::sq;
using especia::sqrt_of_pi;


/**
 * The Gaussian.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the Gaussian at @c x.
 */
static R_type f_g(const R_type &x, const R_type &gamma) {
    return (1.0 / (sqrt_of_pi * gamma)) * exp(-sq(x / gamma));
}

/**
 * The Lorentzian.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the Lorentzian at @c x.
 */
static R_type f_l(const R_type &x, const R_type &gamma) {
    return 1.0 / ((pi * gamma) * (1.0 + sq(x / gamma)));
}

/**
 * The irrational function used in the extended pseudo-Voigt approximation.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the function at @c x.
 */
static R_type f_i(const R_type &x, const R_type &gamma) {
    return 1.0 / ((2.0 * gamma) * pow(1.0 + sq(x / gamma), 1.5));
}

/**
 * The squared hyperbolic secant function used in the extended pseudo-Voigt approximation.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the function at @c x.
 */
static R_type f_p(const R_type &x, const R_type &gamma) {
    return 1.0 / (2.0 * gamma * sq(cosh(x / gamma)));
}

/**
 * An univariate polynomial of degree 6.
 *
 * @param[in] x The abscissa value.
 * @param[in] h0 The coefficient for the monomial of degree 0.
 * @param[in] h1 The coefficient for the monomial of degree 1.
 * @param[in] h2 The coefficient for the monomial of degree 2.
 * @param[in] h3 The coefficient for the monomial of degree 3.
 * @param[in] h4 The coefficient for the monomial of degree 4.
 * @param[in] h5 The coefficient for the monomial of degree 5.
 * @param[in] h6 The coefficient for the monomial of degree 6.
 * @return the value of the polynomial at @c x.
 */
static R_type poly(const R_type &x,
            const R_type &h0,
            const R_type &h1,
            const R_type &h2,
            const R_type &h3,
            const R_type &h4,
            const R_type &h5,
            const R_type &h6) {
    return h0 + x * (h1 + x * (h2 + x * (h3 + x * (h4 + x * (h5 + x * h6)))));
}

static R_type poly_w_g(const R_type &r) {
    return 1.0 - r * poly(r, 0.66000, 0.15021, -1.24984, 4.74052, -9.48291, 8.48252, -2.95553);
}

static R_type poly_w_l(const R_type &r) {
    return 1.0 - (1.0 - r) * poly(r, -0.42179, -1.25693, 10.30003, -23.45651, 29.14158, -16.50453, 3.19974);
}

static R_type poly_w_i(const R_type &r) {
    return poly(r, 1.19913, 1.43021, -15.36331, 47.06071, -73.61822, 57.92559, -17.80614);
}

static R_type poly_w_p(const R_type &r) {
    return poly(r, 1.10186, -0.47745, -0.68688, 2.76622, -4.55466, 4.05475, -1.26571);
}

static R_type poly_eta_l(const R_type &r) {
    return r * (1.0 + (1.0 - r) * poly(r, -0.30165, -1.38927, 9.31550, -24.10743, 34.96491, -21.18862, 3.70290));
}

static R_type poly_eta_i(const R_type &r) {
    return (r * (1.0 - r)) * poly(r, 0.25437, -0.14107, 3.23653, -11.09215, 22.10544, -24.12407, 9.76947);
}

static R_type poly_eta_p(const R_type &r) {
    return (r * (1.0 - r)) * poly(r, 1.01579, 1.50429, -9.21815, 23.59717, -39.71134, 32.83023, -10.02142);
}


especia::Pseudo_Voigt::Pseudo_Voigt(const R_type &b, const R_type &d)
        : u((c_g * b) / (c_l * d)),
          r(1.0 / pow(1.0 + u * (0.07842 + u * (4.47163 + u * (2.42843 + u * (u + 2.69269)))), 0.2)),
          gamma_g((c_l * d) / (c_g * r)),
          gamma_l((c_l * d) / (c_l * r)),
          eta(r * (1.36603 - r * (0.47719 - r * 0.11116))) {
}

especia::Pseudo_Voigt::~Pseudo_Voigt() {
}

R_type especia::Pseudo_Voigt::operator()(const R_type &x) const {
    return (1.0 - eta) * f_g(x, gamma_g) + eta * f_l(x, gamma_l);
}

const R_type especia::Pseudo_Voigt::c_g = 2.0 * sqrt(log(2.0));
const R_type especia::Pseudo_Voigt::c_l = 2.0;


especia::Extended_Pseudo_Voigt::Extended_Pseudo_Voigt(const R_type &b, const R_type &d)
        : u(c_g * b + c_l * d),
          r(c_l * d / u),
          gamma_g(u * poly_w_g(r) / c_g),
          gamma_l(u * poly_w_l(r) / c_l),
          gamma_i(u * poly_w_i(r) / c_i),
          gamma_p(u * poly_w_p(r) / c_p),
          eta_l(poly_eta_l(r)),
          eta_i(poly_eta_i(r)),
          eta_p(poly_eta_p(r)) {
}

especia::Extended_Pseudo_Voigt::~Extended_Pseudo_Voigt() {
}

R_type especia::Extended_Pseudo_Voigt::operator()(const R_type &x) const {
    return (1.0 - eta_l - eta_i - eta_p) * f_g(x, gamma_g) +
           eta_l * f_l(x, gamma_l) +
           eta_i * f_i(x, gamma_i) +
           eta_p * f_p(x, gamma_p);
}

const R_type especia::Extended_Pseudo_Voigt::c_g = 2.0 * sqrt(log(2.0));
const R_type especia::Extended_Pseudo_Voigt::c_l = 2.0;
const R_type especia::Extended_Pseudo_Voigt::c_i = 2.0 * sqrt(pow(2.0, 2.0 / 3.0) - 1.0);
const R_type especia::Extended_Pseudo_Voigt::c_p = 2.0 * log(sqrt(2.0) + 1.0);


especia::Many_Multiplet::Many_Multiplet()
        : u(0.0), z(1.0), c(0.0), b(1.0), a(0.0) {
}

especia::Many_Multiplet::Many_Multiplet(const R_type q[])
        : u(1.0E+08 / (1.0E+08 / q[0] + q[6] * (q[7] * micro) * (q[7] * micro + 2.0))),
          z((1.0 + q[2]) * (1.0 + q[3] / c0)),
          c(u * z),
          b(q[4] * c / c0),
          a(c1 * q[1] * pow(10.0, q[5]) * (u * c)) {
}

especia::Many_Multiplet::~Many_Multiplet() {
}

R_type especia::Many_Multiplet::operator()(const R_type &x) const {
    return a * truncate(f_g, x - c, b, 4.0);
}

const R_type especia::Many_Multiplet::c0 = 1.0E-03 * speed_of_light;
const R_type especia::Many_Multiplet::c1 = 1.0E-06 * sq(elementary_charge) /
                                           (4.0 * electric_constant * electron_mass * sq(speed_of_light));


especia::Intergalactic_Doppler::Intergalactic_Doppler()
        : z(1.0), c(0.0), b(1.0), a(0.0) {
}

especia::Intergalactic_Doppler::Intergalactic_Doppler(const R_type q[])
        : z((1.0 + q[2]) * (1.0 + q[3] / c0)),
          c(q[0] * z),
          b(q[4] * c / c0),
          a(c1 * q[1] * pow(10.0, q[5]) * (q[0] * c)) {
}

especia::Intergalactic_Doppler::~Intergalactic_Doppler() {
}

R_type especia::Intergalactic_Doppler::operator()(const R_type &x) const {
    return a * truncate(f_g, x - c, b, 4.0);
}

const R_type especia::Intergalactic_Doppler::c0 = 1.0E-03 * speed_of_light;
const R_type especia::Intergalactic_Doppler::c1 = 1.0E-06 * sq(elementary_charge) /
                                                  (4.0 * electric_constant * electron_mass * sq(speed_of_light));
