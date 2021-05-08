//! @file decompose.cxx
//! Symmetric eigenproblem solvers calling the LAPACK routines.

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
#include <algorithm>

#include "decompose.h"

using std::copy;
using std::invalid_argument;
using std::max;
using std::runtime_error;
using std::string;
using std::swap;
using std::valarray;

using especia::real;
using especia::integer;


#define LAPACK_NAME_DOUBLE(x) d##x##_
#define LAPACK_NAME_SINGLE(x) s##x##_
#define LAPACK_NAME_R_TYPE(x) LAPACK_NAME_DOUBLE(x)

extern "C" {
/// Interface to LAPACK routine @c [DS]LAMCH.
real LAPACK_NAME_R_TYPE(lamch)(const char &cmach);

/// Interface to LAPACK routine @c [DS]SYEVD.
void LAPACK_NAME_R_TYPE(syevd)(const char &job,
                               const char &uplo,
                               const integer &n,
                               real A[],
                               const integer &lda,
                               real w[],
                               real work[],
                               const integer &lwork,
                               integer iwork[],
                               const integer &liwork,
                               integer &info);

/// Interface to LAPACK routine @c [DS]SYEVR.
void LAPACK_NAME_R_TYPE(syevr)(const char &job,
                               const char &range,
                               const char &uplo,
                               const integer &n,
                               real A[],
                               const integer &lda,
                               const real &vl,
                               const real &vu,
                               const integer &il,
                               const integer &iu,
                               const real &abstol,
                               integer &m,
                               real w[],
                               real Z[],
                               const integer &ldz,
                               integer isupp[],
                               real work[],
                               const integer &lwork,
                               integer iwork[],
                               const integer &liwork,
                               integer &info);

/// Interface to LAPACK routine @c [DS]SYEVX.
void LAPACK_NAME_R_TYPE(syevx)(const char &job,
                               const char &range,
                               const char &uplo,
                               const integer &n,
                               real A[],
                               const integer &lda,
                               const real &vl,
                               const real &vu,
                               const integer &il,
                               const integer &iu,
                               const real &abstol,
                               integer &m,
                               real w[],
                               real Z[],
                               const integer &ldz,
                               real work[],
                               const integer &lwork,
                               integer iwork[],
                               integer ifail[],
                               integer &info);
}

/// The LAPACK job parameter (here: compute eigenvalues and eigenvectors).
static const char job = 'V';

/// The LAPACK range parameter (here: compute all eigenvalues and eigenvectors).
static const char range = 'A';

/// The LAPACK matrix store parameter (here: use the upper triangular part)
static const char uplo = 'U';

/// The LAPACK lower range limit (here: not used).
static const real vl = 0.0;

/// The LAPACK upper range limit (here: not used).
static const real vu = 0.0;

/// The LAPACK lower range index (here: not used).
static const integer il = 0;

/// The LAPACK upper range index (here: not used).
static const integer iu = 0;


especia::D_Decompose::D_Decompose(natural m) : n(integer(m)), work(), iwork() {
    lapack_inquire(n, lwork, liwork);

    work.resize(static_cast<size_t>(lwork));
    iwork.resize(static_cast<size_t>(liwork));
}

especia::D_Decompose::~D_Decompose() = default;

void especia::D_Decompose::operator()(const real A[], real Z[], real w[]) const {
    copy(&A[0], &A[n * n], Z);

    lapack_do(Z, w);
}

void especia::D_Decompose::lapack_do(real Z[], real w[]) const {
    integer info = 0;

    LAPACK_NAME_R_TYPE(syevd)(job, uplo, n, &Z[0], n, w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // ok
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::D_Decompose::lapack_inquire(integer n, integer &lwork, integer &liwork) {
    integer info;
    real work;

    LAPACK_NAME_R_TYPE(syevd)(job, uplo, n, nullptr, n, nullptr, &work, -1, &liwork, -1, info);

    if (info == 0) {
        lwork = static_cast<integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const string especia::D_Decompose::message_int_err = // NOLINT
        "especia::D_Decompose() Error: internal error in LAPACK";
const string especia::D_Decompose::message_ill_arg = // NOLINT
        "especia::D_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::R_Decompose::R_Decompose(natural m)
        : n(integer(m)), work(), iwork(), isupp(2 * max<natural>(1, m)), awork(m * m) {
    lapack_inquire(n, lwork, liwork);

    work.resize(static_cast<size_t>(lwork));
    iwork.resize(static_cast<size_t>(liwork));
}

especia::R_Decompose::~R_Decompose() = default;

void especia::R_Decompose::operator()(const real A[], real Z[], real w[]) const {
    copy(&A[0], &A[n * n], &awork[0]);

    lapack_do(Z, w);
}

void especia::R_Decompose::lapack_do(real Z[], real w[]) const {
    integer m = 0;
    integer info = 0;

    LAPACK_NAME_R_TYPE(syevr)(job, range, uplo, n, &awork[0], n, vl, vu, il, iu, abstol, m, w, Z, n,
                              &isupp[0], &work[0], lwork, &iwork[0], liwork,
                              info);

    if (info == 0) {
        // ok
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::R_Decompose::lapack_inquire(integer n, integer &lwork, integer &liwork) {
    integer info;
    integer m;
    real work;

    LAPACK_NAME_R_TYPE(syevr)(job, range, uplo, n, nullptr, n, vl, vu, il, iu, abstol, m, nullptr, nullptr, n,
                              nullptr, &work, -1, &liwork, -1,
                              info);

    if (info == 0) {
        lwork = static_cast<integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const real especia::R_Decompose::abstol = LAPACK_NAME_R_TYPE(lamch)('S'); // NOLINT
const string especia::R_Decompose::message_int_err = // NOLINT
        "especia::R_Decompose() Error: internal error in LAPACK";
const string especia::R_Decompose::message_ill_arg = // NOLINT
        "especia::R_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::X_Decompose::X_Decompose(natural m) : n(integer(m)), work(), iwork(5 * m), ifail(m), awork(m * m) {
    lapack_inquire(n, lwork);

    work.resize(static_cast<size_t>(lwork));
}

especia::X_Decompose::~X_Decompose() = default;

void especia::X_Decompose::operator()(const real A[], real Z[], real w[]) const {
    copy(&A[0], &A[n * n], &awork[0]);

    lapack_do(Z, w);
}

void especia::X_Decompose::lapack_do(real Z[], real w[]) const {
    integer m = 0;
    integer info = 0;

    LAPACK_NAME_R_TYPE(syevx)(job, range, uplo, n, &awork[0], n, vl, vu, il, iu, abstol, m, w, Z, n,
                              &work[0], lwork, &iwork[0], &ifail[0],
                              info);

    if (info == 0) {
        // ok
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::X_Decompose::lapack_inquire(integer n, integer &lwork) {
    integer info;
    integer m;
    real work;

    LAPACK_NAME_R_TYPE(syevx)(job, range, uplo, n, nullptr, n, vl, vu, il, iu, abstol, m, nullptr, nullptr, n,
                              &work, -1, nullptr, nullptr,
                              info);

    if (info == 0) {
        lwork = static_cast<integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const real especia::X_Decompose::abstol = real(2) * LAPACK_NAME_R_TYPE(lamch)('S'); // NOLINT
const string especia::X_Decompose::message_int_err = // NOLINT
        "especia::X_Decompose() Error: internal error in LAPACK";
const string especia::X_Decompose::message_ill_arg = // NOLINT
        "especia::X_Decompose() Error: illegal argument(s) in call to LAPACK";

#undef LAPACK_NAME_R_TYPE
#undef LAPACK_NAME_SINGLE
#undef LAPACK_NAME_DOUBLE
