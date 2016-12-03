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
using std::pow;

using especia::PI;
using especia::sqr;
using especia::SQRT_OF_PI;

inline
double f_g(const double &x, const double &gamma) {
    return (1.0 / (SQRT_OF_PI * gamma)) * exp(-sqr(x / gamma));
}

inline
double f_l(const double &x, const double &gamma) {
    return 1.0 / ((PI * gamma) * (1.0 + sqr(x / gamma)));
}

inline
double f_i(const double &x, const double &gamma) {
    return 1.0 / ((2.0 * gamma) * pow(1.0 + sqr(x / gamma), 1.5));
}

inline
double f_p(const double &x, const double &gamma) {
    return 1.0 / (2.0 * gamma * sqr(cosh(x / gamma)));
}

inline
double poly(const double &rho,
            const double &h0,
            const double &h1,
            const double &h2,
            const double &h3,
            const double &h4,
            const double &h5,
            const double &h6) {
    return h0 + rho * (h1 + rho * (h2 + rho * (h3 + rho * (h4 + rho * (h5 + rho * h6)))));
}

inline
double w_g(const double &rho) {
    return 1.0 - rho * poly(rho, 0.66000, 0.15021, -1.24984, 4.74052, -9.48291, 8.48252, -2.95553);
}

inline
double w_l(const double &rho) {
    return 1.0 - (1.0 - rho) * poly(rho, -0.42179, -1.25693, 10.30003, -23.45651, 29.14158, -16.50453, 3.19974);
}

inline
double w_i(const double &rho) {
    return poly(rho, 1.19913, 1.43021, -15.36331, 47.06071, -73.61822, 57.92559, -17.80614);
}

inline
double w_p(const double &rho) {
    return poly(rho, 1.10186, -0.47745, -0.68688, 2.76622, -4.55466, 4.05475, -1.26571);
}

inline
double eta_l(const double &rho) {
    return rho * (1.0 + (1.0 - rho) * poly(rho, -0.30165, -1.38927, 9.31550, -24.10743, 34.96491, -21.18862, 3.70290));
}

inline
double eta_i(const double &rho) {
    return (rho * (1.0 - rho)) * poly(rho, 0.25437, -0.14107, 3.23653, -11.09215, 22.10544, -24.12407, 9.76947);
}

inline
double eta_p(const double &rho) {
    return (rho * (1.0 - rho)) * poly(rho, 1.01579, 1.50429, -9.21815, 23.59717, -39.71134, 32.83023, -10.02142);
}

double pseudo_voigt(const double &x, const double &gamma_g, const double &gamma_l, const double &eta) {
    return (1.0 - eta) * f_g(x, gamma_g) + eta * f_l(x, gamma_l);
}

double extended_pseudo_voigt(const double &x,
                             const double &gamma_g,
                             const double &gamma_l,
                             const double &gamma_i,
                             const double &gamma_p,
                             const double &eta_l,
                             const double &eta_i,
                             const double &eta_p) {
    return (1.0 - eta_l - eta_i - eta_p) * f_g(x, gamma_g) + eta_l * f_l(x, gamma_l) + eta_i * f_i(x, gamma_i) +
           eta_p * f_p(x, gamma_p);
}

especia::doppler_mm::doppler_mm() : b(1.0), c(0.0), y(0.0) {
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
    const double x = 1.0e+08 / (1.0e+08 / a[0] + a[6] * v);

    y = x * (1.0 + a[2]) * (1.0 + a[3] / SPEED_OF_LIGHT);
    b = a[4] * y / SPEED_OF_LIGHT;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (x * y);
}

especia::doppler_ig::doppler_ig() : b(1.0), c(0.0), y(0.0) {
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
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / SPEED_OF_LIGHT);
    b = a[4] * y / SPEED_OF_LIGHT;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
}

especia::doppler_is::doppler_is() : b(1.0), c(0.0), y(0.0) {
}

especia::doppler_is::doppler_is(const double a[]) {
    assign(a);
}

especia::doppler_is::~doppler_is() {
}

double especia::doppler_is::operator()(double x) const {
    using std::abs;

    return (abs(x - y) < 4.0 * b) ? c * f_g(x - y, b) : 0.0;
}

void especia::doppler_is::assign(const double a[]) {
    y = a[0] * (1.0 + a[2] / SPEED_OF_LIGHT);
    b = a[3] * y / SPEED_OF_LIGHT;
    c = 8.85280e-21 * a[1] * pow(10.0, a[4]) * a[0];
}

especia::voigt_ig::voigt_ig() : b(1.0), c(0.0), d(1.0), y(0.0), z(0.0) {
}

especia::voigt_ig::voigt_ig(const double a[]) {
    assign(a);
}

especia::voigt_ig::~voigt_ig() {
}

double especia::voigt_ig::operator()(double x) const {
    return c * pseudo_voigt(x - y, b, d, z);
}

void especia::voigt_ig::assign(const double a[]) {
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / SPEED_OF_LIGHT);
    b = a[4] * y / SPEED_OF_LIGHT;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
    d = 2.65442e-20 * a[6] * (a[0] * y);

    z = SQRT_OF_LN_OF_2 * (b / d);
    z = 1.0 / pow(1.0 + z * (0.07842 + z * (4.47163 + z * (2.42843 + z * (z + 2.69269)))), 0.2);

    d = d / z;
    b = d / SQRT_OF_LN_OF_2;

    z = z * (1.36603 - z * (0.47719 - z * 0.11116));
}

// References
//
// T. Ida, M. Ando, H. Toraya (2000)
//   Extended pseudo-Voigt function for approximating the Voigt profile,
//   J. Appl. Chryst., 33, 1311, ISSN 0021-8898
