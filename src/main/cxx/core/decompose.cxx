/// @file decompose.cxx
/// Symmetric eigenproblem solvers calling the LAPACK routines.
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
#include <algorithm>

#include "decompose.h"

#define LAPACK_NAME(x) d##x##_

using std::copy;
using std::max;
using std::runtime_error;
using std::string;
using std::swap;
using std::valarray;

using especia::R_type;
using especia::Z_type;

/**
 * Interface to LAPACK eigenvalue routines (version 3.0).
 */
extern "C" {
R_type LAPACK_NAME(lamch)(const char &cmach);

void LAPACK_NAME(syevd)(const char &job,
                        const char &uplo,
                        const Z_type &n,
                        R_type A[], const Z_type &lda,
                        R_type w[],
                        R_type work[], const Z_type &lwork,
                        Z_type iwork[], const Z_type &liwork,
                        Z_type &info);

void LAPACK_NAME(syevr)(const char &job,
                        const char &range,
                        const char &uplo,
                        const Z_type &n,
                        R_type A[], const Z_type &lda,
                        const R_type &vl, const R_type &vu,
                        const Z_type &il, const Z_type &iu,
                        const R_type &abstol,
                        Z_type &m,
                        R_type w[],
                        R_type Z[], const Z_type &ldz,
                        Z_type isupp[],
                        R_type work[], const Z_type &lwork,
                        Z_type iwork[], const Z_type &liwork,
                        Z_type &info);

void LAPACK_NAME(syevx)(const char &job,
                        const char &range,
                        const char &uplo,
                        const Z_type &n,
                        R_type A[], const Z_type &lda,
                        const R_type &vl, const R_type &vu,
                        const Z_type &il, const Z_type &iu,
                        const R_type &abstol,
                        Z_type &m,
                        R_type w[],
                        R_type Z[], const Z_type &ldz,
                        R_type work[], const Z_type &lwork,
                        Z_type iwork[],
                        Z_type ifail[], Z_type &info);
}

const R_type safe_minimum = LAPACK_NAME(lamch)('s');

const string especia::D_Decompose::message_int_err = "especia::D_Decompose() Error: internal error in LAPACK";
const string especia::D_Decompose::message_ill_arg = "especia::D_Decompose() Error: illegal argument(s) in call to LAPACK";

especia::D_Decompose::D_Decompose(N_type n)
        : job('V'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::D_Decompose::~D_Decompose() {
}

void especia::D_Decompose::operator()(N_type n_in, const R_type A[], R_type Z[], R_type w[]) throw(runtime_error) {
    copy(&A[0], &A[n_in * n_in], Z);

    if (n_in != static_cast<N_type>(n)) {
        resize_workspace(n_in);
    }

    // The regular call.
    LAPACK_NAME(syevd)(job, uplo, n, &Z[0], max(1, n), w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::D_Decompose::resize_workspace(N_type n_in) {
    n = static_cast<Z_type>(n_in);

    // The workspace query.
    LAPACK_NAME(syevd)(job, uplo, n, 0, max(1, n), 0, &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        liwork = iwork[0];
        work.resize(static_cast<N_type>(lwork));
        iwork.resize(static_cast<N_type>(liwork));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::D_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}


const string especia::R_Decompose::message_int_err = "especia::R_Decompose() Error: internal error in LAPACK";
const string especia::R_Decompose::message_ill_arg = "especia::R_Decompose() Error: illegal argument(s) in call to LAPACK";

especia::R_Decompose::R_Decompose(N_type n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::R_Decompose::~R_Decompose() {
}

void especia::R_Decompose::operator()(N_type n_in, const R_type A[], R_type Z[], R_type w[]) throw(runtime_error) {
    valarray<R_type> C(A, n_in * n_in);

    if (n_in != static_cast<N_type>(n)) {
        resize_workspace(n_in);
    }

    // The regular call.
    LAPACK_NAME(syevr)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, safe_minimum, m, w, Z,
                       max(1, n), &isupp[0], &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::R_Decompose::resize_workspace(N_type n_in) {
    n = static_cast<Z_type>(n_in);

    // The workspace query.
    LAPACK_NAME(syevr)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, safe_minimum, m, 0, 0,
                       max(1, n), &isupp[0], &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        liwork = iwork[0];
        work.resize(static_cast<N_type>(lwork));
        iwork.resize(static_cast<N_type>(liwork));
        isupp.resize(static_cast<N_type>(2 * max(1, n)));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::R_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}


const string especia::X_Decompose::message_int_err = "especia::X_Decompose() Error: internal error in LAPACK";
const string especia::X_Decompose::message_ill_arg = "especia::X_Decompose() Error: illegal argument(s) in call to LAPACK";

especia::X_Decompose::X_Decompose(N_type n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(), ifail() {
    resize_workspace(n);
}

especia::X_Decompose::~X_Decompose() {
}

void especia::X_Decompose::operator()(N_type n_in, const R_type A[], R_type Z[], R_type w[]) throw(runtime_error) {
    valarray<R_type> C(A, n_in * n_in);

    if (n_in != static_cast<N_type>(n)) {
        resize_workspace(n_in);
    }

    // The regular call.
    LAPACK_NAME(syevx)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum, m, w, Z,
                       max(1, n), &work[0], lwork, &iwork[0], &ifail[0], info);

    if (info == 0) {
        // Transform from column-major into row-major layout.
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::X_Decompose::resize_workspace(N_type n_in) {
    n = static_cast<Z_type>(n_in);

    // The workspace query.
    LAPACK_NAME(syevx)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum, m, 0, 0,
                       max(1, n), &work[0], -1, &iwork[0], &ifail[0], info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        work.resize(static_cast<N_type>(lwork));
        iwork.resize(static_cast<N_type>(5 * n));
        ifail.resize(static_cast<N_type>(n));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw runtime_error(message_ill_arg);
    }
}

void especia::X_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < n; ++i, i0 += n) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n) {
            swap(A[ij], A[ji]);
        }
    }
}

#undef LAPACK_NAME