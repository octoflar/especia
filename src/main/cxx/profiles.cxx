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
using std::exp;
using std::pow;

const double speed_of_light = 299792.458;

double
RQ::voigt(double x, double b, double d) {
    const double a = sqrt_of_ln_of_2 * (b / d);
    const double c = pow(1.0 + a * (0.07842 + a * (4.47163 + a * (2.42843 +
                                                                  a * (a + 2.69269)))), -0.2);
    const double z = a * (1.36603 - a * (0.47719 - a * 0.11116));

    d = d / c;
    b = d / sqrt_of_ln_of_2;

    return pseudovoigt(x, b, d, z);
}

RQ::doppler_mm_pf::doppler_mm_pf() : b(1.0), c(0.0), y(0.0) {
}

RQ::doppler_mm_pf::doppler_mm_pf(const double a[]) {
    assign(a);
}

RQ::doppler_mm_pf::~doppler_mm_pf() {
}

void
RQ::doppler_mm_pf::assign(const double a[]) {
    const double u = 1.0e-05 * a[7];
    const double v = u * (u + 2.0);
    const double x = 1.0e+08 / (1.0e+08 / a[0] + a[6] * v);

    y = x * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (x * y);
}

RQ::doppler_pf::doppler_pf() : b(1.0), c(0.0), y(0.0) {
}

RQ::doppler_pf::doppler_pf(const double a[]) {
    assign(a);
}

RQ::doppler_pf::~doppler_pf() {
}

void
RQ::doppler_pf::assign(const double a[]) {
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
}

RQ::doppler_is_pf::doppler_is_pf() : b(1.0), c(0.0), y(0.0) {
}

RQ::doppler_is_pf::doppler_is_pf(const double a[]) {
    assign(a);
}

RQ::doppler_is_pf::~doppler_is_pf() {
}

void
RQ::doppler_is_pf::assign(const double a[]) {
    y = a[0] * (1.0 + a[2] / speed_of_light);
    b = a[3] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[4]) * a[0];
}

RQ::voigt_pf::voigt_pf() : b(1.0), c(0.0), d(1.0), y(0.0), z(0.0) {
}

RQ::voigt_pf::voigt_pf(const double a[]) {
    assign(a);
}

RQ::voigt_pf::~voigt_pf() {
}

void
RQ::voigt_pf::assign(const double a[]) {
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
    d = 2.65442e-20 * a[6] * (a[0] * y);

    z = sqrt_of_ln_of_2 * (b / d);
    z = pow(1.0 + z * (0.07842 + z * (4.47163 + z * (2.42843 + z * (z + 2.69269)))), -0.2);

    d = d / z;
    b = d / sqrt_of_ln_of_2;

    z = z * (1.36603 - z * (0.47719 - z * 0.11116));
}

// References
//
// T. Ida, M. Ando, H. Toraya (2000)
//   Extended pseudo-Voigt function for approximating the Voigt profile,
//   J. Appl. Chryst., 33, 1311, ISSN 0021-8898
