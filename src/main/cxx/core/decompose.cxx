/// @file decompose.cxx
/// Symmetric eigenproblem solvers calling the LAPACK routines.
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
#include <algorithm>

#include "decompose.h"

#define LAPACK_NAME(x) d##x##_

using std::copy;
using std::max;
using std::runtime_error;
using std::string;
using std::swap;
using std::valarray;

using especia::Real_t;
using especia::Zint_t;

/**
 * Interface to LAPACK eigenvalue routines (version 3.0).
 */
extern "C" {
Real_t LAPACK_NAME(lamch)(const char &cmach);

void LAPACK_NAME(syevd)(const char &job,
                        const char &uplo,
                        const Zint_t &n,
                        Real_t A[], const Zint_t &lda,
                        Real_t w[],
                        Real_t work[], const Zint_t &lwork,
                        Zint_t iwork[], const Zint_t &liwork,
                        Zint_t &info);

void LAPACK_NAME(syevr)(const char &job,
                        const char &range,
                        const char &uplo,
                        const Zint_t &n,
                        Real_t A[], const Zint_t &lda,
                        const Real_t &vl, const Real_t &vu,
                        const Zint_t &il, const Zint_t &iu,
                        const Real_t &abstol,
                        Zint_t &m,
                        Real_t w[],
                        Real_t Z[], const Zint_t &ldz,
                        Zint_t isupp[],
                        Real_t work[], const Zint_t &lwork,
                        Zint_t iwork[], const Zint_t &liwork,
                        Zint_t &info);

void LAPACK_NAME(syevx)(const char &job,
                        const char &range,
                        const char &uplo,
                        const Zint_t &n,
                        Real_t A[], const Zint_t &lda,
                        const Real_t &vl, const Real_t &vu,
                        const Zint_t &il, const Zint_t &iu,
                        const Real_t &abstol,
                        Zint_t &m,
                        Real_t w[],
                        Real_t Z[], const Zint_t &ldz,
                        Real_t work[], const Zint_t &lwork,
                        Zint_t iwork[],
                        Zint_t ifail[], Zint_t &info);
}

const Real_t safe_minimum = LAPACK_NAME(lamch)('s');

const string especia::D_Decompose::int_err = "especia::D_Decompose(): Error: internal error in LAPACK routine DSYEVD";
const string especia::D_Decompose::ill_arg = "especia::D_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVD";

especia::D_Decompose::D_Decompose(Nint_t n)
        : job('V'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::D_Decompose::~D_Decompose() {
}

void especia::D_Decompose::operator()(Nint_t k, const Real_t A[], Real_t Z[], Real_t w[]) throw(runtime_error) {
    copy(&A[0], &A[k * k], Z);

    if (k != static_cast<Nint_t>(n)) {
        resize_workspace(k);
    }

    // The regular call.
    LAPACK_NAME(syevd)(job, uplo, n, &Z[0], max(1, n), w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::D_Decompose::resize_workspace(Nint_t k) {
    n = static_cast<Zint_t>(k);

    // The workspace query.
    LAPACK_NAME(syevd)(job, uplo, n, 0, max(1, n), 0, &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<Zint_t>(work[0]);
        liwork = iwork[0];
        work.resize(static_cast<Nint_t>(lwork));
        iwork.resize(static_cast<Nint_t>(liwork));
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::D_Decompose::transpose(Real_t A[]) const {
    for (Zint_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Zint_t j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}


const string especia::R_Decompose::int_err = "especia::R_Decompose(): Error: internal error in LAPACK routine DSYEVR";
const string especia::R_Decompose::ill_arg = "especia::R_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVR";

especia::R_Decompose::R_Decompose(Nint_t n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::R_Decompose::~R_Decompose() {
}

void especia::R_Decompose::operator()(Nint_t k, const Real_t A[], Real_t Z[], Real_t w[]) throw(runtime_error) {
    valarray<Real_t> C(A, k * k);

    if (k != static_cast<Nint_t>(n)) {
        resize_workspace(k);
    }

    // The regular call.
    LAPACK_NAME(syevr)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, safe_minimum, m, w, Z,
                       max(1, n), &isupp[0], &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::R_Decompose::resize_workspace(Nint_t k) {
    n = static_cast<Zint_t>(k);

    // The workspace query.
    LAPACK_NAME(syevr)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, safe_minimum, m, 0, 0,
                       max(1, n), &isupp[0], &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<Zint_t>(work[0]);
        liwork = iwork[0];
        work.resize(static_cast<Nint_t>(lwork));
        iwork.resize(static_cast<Nint_t>(liwork));
        isupp.resize(static_cast<Nint_t>(2 * max(1, n)));
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::R_Decompose::transpose(Real_t A[]) const {
    for (Zint_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Zint_t j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}


const string especia::X_Decompose::int_err = "especia::X_Decompose(): Error: internal error in LAPACK routine DSYEVX";
const string especia::X_Decompose::ill_arg = "especia::X_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVX";

especia::X_Decompose::X_Decompose(Nint_t n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(), ifail() {
    resize_workspace(n);
}

especia::X_Decompose::~X_Decompose() {
}

void especia::X_Decompose::operator()(Nint_t k, const Real_t A[], Real_t Z[], Real_t w[]) throw(runtime_error) {
    valarray<Real_t> C(A, k * k);

    if (k != static_cast<Nint_t>(n)) {
        resize_workspace(k);
    }

    // The regular call.
    LAPACK_NAME(syevx)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum, m, w, Z,
                       max(1, n), &work[0], lwork, &iwork[0], &ifail[0], info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::X_Decompose::resize_workspace(Nint_t k) {
    n = static_cast<Zint_t>(k);

    // The workspace query.
    LAPACK_NAME(syevx)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum, m, 0, 0,
                       max(1, n), &work[0], -1, &iwork[0], &ifail[0], info);

    if (info == 0) {
        lwork = static_cast<Zint_t>(work[0]);
        work.resize(static_cast<Nint_t>(lwork));
        iwork.resize(static_cast<Nint_t>(5 * n));
        ifail.resize(static_cast<Nint_t>(n));
    } else if (info > 0) {
        throw runtime_error(int_err);
    } else {
        throw runtime_error(ill_arg);
    }
}

void especia::X_Decompose::transpose(Real_t A[]) const {
    for (Zint_t i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Zint_t j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}

#undef LAPACK_NAME