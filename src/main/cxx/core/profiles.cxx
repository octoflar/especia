//! @file profiles.cxx
//! Profile functions.

// Copyright (c) 2021. Ralf Quast
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
#include "profiles.h"

using std::abs;
using std::cosh;
using std::exp;
using std::log;
using std::pow;
using std::sqrt;

using especia::real;
using especia::micro;
using especia::pi;
using especia::sq;
using especia::sqrt_of_ln_two;
using especia::sqrt_of_pi;


/**
 * The Gaussian.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the Gaussian at @c x.
 */
static real f_g(const real &x, const real &gamma) {
    return (1.0 / (sqrt_of_pi * gamma)) * exp(-sq(x / gamma));
}

/**
 * The Lorentzian.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the Lorentzian at @c x.
 */
static real f_l(const real &x, const real &gamma) {
    return 1.0 / ((pi * gamma) * (1.0 + sq(x / gamma)));
}

/**
 * The irrational function used in the extended pseudo-Voigt approximation.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the function at @c x.
 */
static real f_i(const real &x, const real &gamma) {
    return 1.0 / ((2.0 * gamma) * pow(1.0 + sq(x / gamma), 1.5));
}

/**
 * The squared hyperbolic secant function used in the extended pseudo-Voigt approximation.
 *
 * @param[in] x The abscissa value (arbitrary unit).
 * @param[in] gamma The width (arbitrary unit).
 * @return the value of the function at @c x.
 */
static real f_p(const real &x, const real &gamma) {
    return 1.0 / (2.0 * gamma * sq(cosh(x / gamma)));
}

template<class F>
static real truncate(const F &f, const real &x, const real &b, const real &c) {
    return abs(x) < c * b ? f(x, b) : real(0.0);
}

template<class T>
static T poly(const T &x, const T &h0, const T &h1, const T &h2, const T &h3, const T &h4, const T &h5, const T &h6) {
    return h0 + x * (h1 + x * (h2 + x * (h3 + x * (h4 + x * (h5 + x * h6)))));
}

static real poly_w_g(const real &r) {
    return 1.0 - r * poly(r, 0.66000, 0.15021, -1.24984, 4.74052, -9.48291, 8.48252, -2.95553);
}

static real poly_w_l(const real &r) {
    return 1.0 - (1.0 - r) * poly(r, -0.42179, -1.25693, 10.30003, -23.45651, 29.14158, -16.50453, 3.19974);
}

static real poly_w_i(const real &r) {
    return poly(r, 1.19913, 1.43021, -15.36331, 47.06071, -73.61822, 57.92559, -17.80614);
}

static real poly_w_p(const real &r) {
    return poly(r, 1.10186, -0.47745, -0.68688, 2.76622, -4.55466, 4.05475, -1.26571);
}

static real poly_eta_l(const real &r) {
    return r * (1.0 + (1.0 - r) * poly(r, -0.30165, -1.38927, 9.31550, -24.10743, 34.96491, -21.18862, 3.70290));
}

static real poly_eta_i(const real &r) {
    return (r * (1.0 - r)) * poly(r, 0.25437, -0.14107, 3.23653, -11.09215, 22.10544, -24.12407, 9.76947);
}

static real poly_eta_p(const real &r) {
    return (r * (1.0 - r)) * poly(r, 1.01579, 1.50429, -9.21815, 23.59717, -39.71134, 32.83023, -10.02142);
}


especia::Pseudo_Voigt::Pseudo_Voigt(const real &b, const real &d)
        : rho(c_g * b / d),
          h(1.0 / pow(1.0 + rho * (0.07842 + rho * (4.47163 + rho * (2.42843 + rho * (rho + 2.69296)))), 0.2)),
          // transposed digits in T. Ida, M. Ando, H. Toraya (2000)                              ^^
          gamma_l(d / h),
          gamma_g(gamma_l / c_g),
          eta(h * (1.36603 - h * (0.47719 - h * 0.11116))) {
}

especia::Pseudo_Voigt::~Pseudo_Voigt() = default;

real especia::Pseudo_Voigt::operator()(const real &x) const {
    return (1.0 - eta) * f_g(x, gamma_g) + eta * f_l(x, gamma_l);
}

const real especia::Pseudo_Voigt::c_g = sqrt_of_ln_two;


especia::Extended_Pseudo_Voigt::Extended_Pseudo_Voigt(const real &b, const real &d)
        : g(c_g * b + d),
          rho(d / g),
          gamma_g(g * poly_w_g(rho) / c_g),
          gamma_l(g * poly_w_l(rho)),
          gamma_i(g * poly_w_i(rho) / c_i),
          gamma_p(g * poly_w_p(rho) / c_p),
          eta_l(poly_eta_l(rho)),
          eta_i(poly_eta_i(rho)),
          eta_p(poly_eta_p(rho)) {
}

especia::Extended_Pseudo_Voigt::~Extended_Pseudo_Voigt() = default;

real especia::Extended_Pseudo_Voigt::operator()(const real &x) const {
    return (1.0 - eta_l - eta_i - eta_p) * f_g(x, gamma_g) +
           eta_l * f_l(x, gamma_l) +
           eta_i * f_i(x, gamma_i) +
           eta_p * f_p(x, gamma_p);
}

const real especia::Extended_Pseudo_Voigt::c_g = sqrt_of_ln_two;
const real especia::Extended_Pseudo_Voigt::c_i = sqrt(pow(2.0, 2.0 / 3.0) - 1.0); // NOLINT
const real especia::Extended_Pseudo_Voigt::c_p = log(sqrt(2.0) + 1.0); // NOLINT


especia::Many_Multiplet::Many_Multiplet()
        : u(0.0), z(1.0), c(0.0), b(0.5), a(1.0) {
}

especia::Many_Multiplet::Many_Multiplet(const real q[])
        : u(1.0E+08 / (1.0E+08 / q[0] + q[6] * (q[7] * micro) * (q[7] * micro + 2.0))),
          z((1.0 + q[2]) * (1.0 + q[3] / c0)),
          c(u * z),
          b(q[4] * c / c0),
          a(c1 * q[1] * pow(10.0, q[5]) * (u * c)) {
}

especia::Many_Multiplet::~Many_Multiplet() = default;

real especia::Many_Multiplet::operator()(const real &x) const {
    return a * truncate(f_g, x - c, b, 4.0);
}

const real especia::Many_Multiplet::c0 = 1.0E-03 * speed_of_light;
const real especia::Many_Multiplet::c1 = 1.0E-06 * sq(elementary_charge) / // NOLINT
                                                   (4.0 * electric_constant * electron_mass * sq(speed_of_light));


especia::Intergalactic_Doppler::Intergalactic_Doppler()
        : z(1.0), c(0.0), b(0.5), a(1.0) {
}

especia::Intergalactic_Doppler::Intergalactic_Doppler(const real q[])
        : z((1.0 + q[2]) * (1.0 + q[3] / c0)),
          c(q[0] * z),
          b(q[4] * c / c0),
          a(c1 * q[1] * pow(10.0, q[5]) * (q[0] * c)) {
}

especia::Intergalactic_Doppler::~Intergalactic_Doppler() = default;

real especia::Intergalactic_Doppler::operator()(const real &x) const {
    return a * truncate(f_g, x - c, b, 4.0);
}

const real especia::Intergalactic_Doppler::c0 = 1.0E-03 * speed_of_light;
const real especia::Intergalactic_Doppler::c1 = 1.0E-06 * sq(elementary_charge) / // NOLINT
                                                          (4.0 * electric_constant * electron_mass * sq(speed_of_light));
