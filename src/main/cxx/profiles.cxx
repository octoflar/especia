// Profile functions
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
#include "profiles.h"

using std::abs;
using std::cosh;
using std::exp;
using std::log;
using std::pow;
using std::sqrt;

using especia::pi;
using especia::sqr;
using especia::sqrt_of_pi;

/**
 * The Gaussian.
 *
 * @param x The abscissa value.
 * @param gamma The width.
 * @return the value of the Gaussian at @code x.
 */
inline
double f_g(const double &x, const double &gamma) {
    return (1.0 / (sqrt_of_pi * gamma)) * exp(-sqr(x / gamma));
}

/**
 * The Lorentzian.
 *
 * @param x The abscissa value.
 * @param gamma The width.
 * @return the value of the Lorentzian at @code x.
 */
inline
double f_l(const double &x, const double &gamma) {
    return 1.0 / ((pi * gamma) * (1.0 + sqr(x / gamma)));
}

/**
 * The irrational function used in the extended pseudo Voigt approximation.
 *
 * @param x The abscissa value.
 * @param gamma The width.
 * @return the value of the function at @code x.
 */
inline
double f_i(const double &x, const double &gamma) {
    return 1.0 / ((2.0 * gamma) * pow(1.0 + sqr(x / gamma), 1.5));
}

/**
 * The squared hyperbolic secant function used in the extended pseudo Voigt approximation.
 *
 * @param x The abscissa value.
 * @param gamma The width.
 * @return the value of the function at @code x.
 */
inline
double f_p(const double &x, const double &gamma) {
    return 1.0 / (2.0 * gamma * sqr(cosh(x / gamma)));
}

/**
 * An univariate polynomial of degree 6.
 *
 * @param x The abscissa value.
 * @param h0 The coefficient for the monomial of degree 0.
 * @param h1 The coefficient for the monomial of degree 1.
 * @param h2 The coefficient for the monomial of degree 2.
 * @param h3 The coefficient for the monomial of degree 3.
 * @param h4 The coefficient for the monomial of degree 4.
 * @param h5 The coefficient for the monomial of degree 5.
 * @param h6 The coefficient for the monomial of degree 6.
 * @return the value of the polynomial at @code rho.
 */
inline
double poly(const double &x,
            const double &h0,
            const double &h1,
            const double &h2,
            const double &h3,
            const double &h4,
            const double &h5,
            const double &h6) {
    return h0 + x * (h1 + x * (h2 + x * (h3 + x * (h4 + x * (h5 + x * h6)))));
}

inline
double poly_w_g(const double &rho) {
    return 1.0 - rho * poly(rho, 0.66000, 0.15021, -1.24984, 4.74052, -9.48291, 8.48252, -2.95553);
}

inline
double poly_w_l(const double &rho) {
    return 1.0 - (1.0 - rho) * poly(rho, -0.42179, -1.25693, 10.30003, -23.45651, 29.14158, -16.50453, 3.19974);
}

inline
double poly_w_i(const double &rho) {
    return poly(rho, 1.19913, 1.43021, -15.36331, 47.06071, -73.61822, 57.92559, -17.80614);
}

inline
double poly_w_p(const double &rho) {
    return poly(rho, 1.10186, -0.47745, -0.68688, 2.76622, -4.55466, 4.05475, -1.26571);
}

inline
double poly_eta_l(const double &rho) {
    return rho * (1.0 + (1.0 - rho) * poly(rho, -0.30165, -1.38927, 9.31550, -24.10743, 34.96491, -21.18862, 3.70290));
}

inline
double poly_eta_i(const double &rho) {
    return (rho * (1.0 - rho)) * poly(rho, 0.25437, -0.14107, 3.23653, -11.09215, 22.10544, -24.12407, 9.76947);
}

inline
double poly_eta_p(const double &rho) {
    return (rho * (1.0 - rho)) * poly(rho, 1.01579, 1.50429, -9.21815, 23.59717, -39.71134, 32.83023, -10.02142);
}


especia::pseudo_voigt::pseudo_voigt(double b, double d) {
    const double p = (c_g * b) / (c_l * d);
    const double q = 1.0 / pow(1.0 + p * (0.07842 + p * (4.47163 + p * (2.42843 + p * (p + 2.69269)))), 0.2);

    gamma_g = (c_l * d) / (q * c_g);
    gamma_l = (c_l * d) / (q * c_l);
    eta = q * (1.36603 - q * (0.47719 - q * 0.11116));
}

especia::pseudo_voigt::~pseudo_voigt() {
}

double especia::pseudo_voigt::operator()(const double &x) const {
    return (1.0 - eta) * f_g(x, gamma_g) + eta * f_l(x, gamma_l);
}

const double especia::pseudo_voigt::c_g = 2.0 * sqrt(log(2.0));
const double especia::pseudo_voigt::c_l = 2.0;


especia::extended_pseudo_voigt::extended_pseudo_voigt(double b, double d) {
    const double sum = c_g * b + c_l * d;
    const double rho = c_l * d / sum;
    gamma_g = sum * poly_w_g(rho) / c_g;
    gamma_l = sum * poly_w_l(rho) / c_l;
    gamma_i = sum * poly_w_i(rho) / c_i;
    gamma_p = sum * poly_w_p(rho) / c_p;
    eta_l = poly_eta_l(rho);
    eta_i = poly_eta_i(rho);
    eta_p = poly_eta_p(rho);
}

especia::extended_pseudo_voigt::~extended_pseudo_voigt() {
}

double especia::extended_pseudo_voigt::operator()(const double &x) const {
    return (1.0 - eta_l - eta_i - eta_p) * f_g(x, gamma_g) +
           eta_l * f_l(x, gamma_l) +
           eta_i * f_i(x, gamma_i) +
           eta_p * f_p(x, gamma_p);
}

const double especia::extended_pseudo_voigt::c_g = 2.0 * sqrt(log(2.0));
const double especia::extended_pseudo_voigt::c_l = 2.0;
const double especia::extended_pseudo_voigt::c_i = 2.0 * sqrt(pow(2.0, 2.0 / 3.0) - 1.0);
const double especia::extended_pseudo_voigt::c_p = 2.0 * log(sqrt(2.0) + 1.0);


especia::doppler_mm::doppler_mm() : y(0.0), b(1.0), c(0.0) {
}

especia::doppler_mm::doppler_mm(const double a[]) {
    assign(a);
}

especia::doppler_mm::~doppler_mm() {
}

double especia::doppler_mm::operator()(double x) const {
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * f_g(x - y, b) : 0.0;
}

void especia::doppler_mm::assign(const double a[]) {
    const double u = 1.0e-05 * a[7];
    const double v = u * (u + 2.0);
    const double w = 1.0e+08 / (1.0e+08 / a[0] + a[6] * v);

    y = w * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (w * y);
}


especia::doppler_ig::doppler_ig() : y(0.0), b(1.0), c(0.0) {
}

especia::doppler_ig::doppler_ig(const double a[]) {
    assign(a);
}

especia::doppler_ig::~doppler_ig() {
}

double especia::doppler_ig::operator()(double x) const {
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * f_g(x - y, b) : 0.0;
}

void especia::doppler_ig::assign(const double a[]) {
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
}
