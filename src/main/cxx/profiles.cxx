// Profile functions
// Copyright (c) 2016, Ralf Quast
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "profiles.h"

using std::abs;
using std::exp;
using std::pow;

const double speed_of_light = 299792.458;

double
RQ::voigt(double x, double b, double d)
{
    const double a = sqrt_of_ln_of_2 * (b / d);
    const double c = pow(1.0 + a * (0.07842 + a * (4.47163 + a * (2.42843 +
        a * (a + 2.69269)))), -0.2);
    const double z = a * (1.36603 - a * (0.47719 - a * 0.11116));

    d = d / c;
    b = d / sqrt_of_ln_of_2;

    return pseudovoigt(x, b, d, z);
}

RQ::gaumm_pf::gaumm_pf() : b(1.0), c(0.0), y(0.0)
{
}

RQ::gaumm_pf::gaumm_pf(const double a[])
{
    assign(a);
}

RQ::gaumm_pf::~gaumm_pf()
{
}

void
RQ::gaumm_pf::assign(const double a[])
{
    const double u = 1.0e-05 * a[7];
    const double v = u * (u + 2.0);
    const double x = 1.0e+08 / (1.0e+08 / a[0] + a[6] * v);

    y = x * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (x * y);
}

RQ::gauss_pf::gauss_pf() : b(1.0), c(0.0), y(0.0)
{
}

RQ::gauss_pf::gauss_pf(const double a[])
{
    assign(a);
}

RQ::gauss_pf::~gauss_pf()
{
}

void
RQ::gauss_pf::assign(const double a[])
{
    y = a[0] * (1.0 + a[2]) * (1.0 + a[3] / speed_of_light);
    b = a[4] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[5]) * (a[0] * y);
}

RQ::test::test() : b(1.0), c(0.0), y(0.0)
{
}

RQ::test::test(const double a[])
{
    assign(a);
}

RQ::test::~test()
{
}

void
RQ::test::assign(const double a[])
{
    y = a[2];
    b = a[3] * y / speed_of_light;
    c = 8.85280e-21 * a[1] * pow(10.0, a[4]) * (a[0] * y);
}

RQ::voigt_pf::voigt_pf() : b(1.0), c(0.0), d(1.0), y(0.0), z(0.0)
{
}

RQ::voigt_pf::voigt_pf(const double a[])
{
    assign(a);
}

RQ::voigt_pf::~voigt_pf()
{
}

void
RQ::voigt_pf::assign(const double a[])
{
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
