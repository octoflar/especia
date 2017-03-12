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

using std::copy;
using std::invalid_argument;
using std::max;
using std::runtime_error;
using std::string;
using std::swap;
using std::valarray;

using especia::R_type;
using especia::Z_type;


#define LAPACK_NAME_DOUBLE(x) d##x##_
#define LAPACK_NAME_SINGLE(x) s##x##_
#define LAPACK_NAME_R_TYPE(x) LAPACK_NAME_DOUBLE(x)

extern "C" {
/**
 * Returns a machine parameter.
 *
 * @param cmach The name of the machine parameter.
 * @return the value of the machine parameter.
 */
R_type LAPACK_NAME_R_TYPE(lamch)(const char &cmach);

/**
 * DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
 * real symmetric matrix A. If eigenvectors are desired, it uses a
 * divide and conquer algorithm.
 *
 * The divide and conquer algorithm makes very mild assumptions about
 * floating point arithmetic. It will work on machines with a guard
 * digit in add/subtract, or on those binary machines without guard
 * digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
 * Cray-2. It could conceivably fail on hexadecimal or decimal machines
 * without guard digits, but we know of none.
 */
void LAPACK_NAME_R_TYPE(syevd)(const char &job,
                               const char &uplo,
                               const Z_type &n,
                               R_type A[],
                               const Z_type &lda,
                               R_type w[],
                               R_type work[],
                               const Z_type &lwork,
                               Z_type iwork[],
                               const Z_type &liwork,
                               Z_type &info);

/**
 * DSYEVR computes selected eigenvalues and, optionally, eigenvectors
 * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 * selected by specifying either a range of values or a range of
 * indices for the desired eigenvalues.
 *
 * Normal execution of DSTEMR may create NaNs and infinities and
 * hence may abort due to a floating point exception in environments
 * which do not handle NaNs and infinities in the ieee standard default
 * manner.
 */
void LAPACK_NAME_R_TYPE(syevr)(const char &job,
                               const char &range,
                               const char &uplo,
                               const Z_type &n,
                               R_type A[],
                               const Z_type &lda,
                               const R_type &vl,
                               const R_type &vu,
                               const Z_type &il,
                               const Z_type &iu,
                               const R_type &abstol,
                               Z_type &m,
                               R_type w[],
                               R_type Z[],
                               const Z_type &ldz,
                               Z_type isupp[],
                               R_type work[],
                               const Z_type &lwork,
                               Z_type iwork[],
                               const Z_type &liwork,
                               Z_type &info);

/**
 * DSYEVX computes selected eigenvalues and, optionally, eigenvectors
 * of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 * selected by specifying either a range of values or a range of indices
 * for the desired eigenvalues.
 */
void LAPACK_NAME_R_TYPE(syevx)(const char &job,
                               const char &range,
                               const char &uplo,
                               const Z_type &n,
                               R_type A[],
                               const Z_type &lda,
                               const R_type &vl,
                               const R_type &vu,
                               const Z_type &il,
                               const Z_type &iu,
                               const R_type &abstol,
                               Z_type &m,
                               R_type w[],
                               R_type Z[],
                               const Z_type &ldz,
                               R_type work[],
                               const Z_type &lwork,
                               Z_type iwork[],
                               Z_type ifail[],
                               Z_type &info);
}

/**
 * The safe minimum, such that its reciprocal does not overflow.
 */
static const R_type safe = LAPACK_NAME_R_TYPE(lamch)('S');


especia::D_Decompose::D_Decompose(N_type n)
        : m(Z_type(n)), work(1), iwork(1) {
    allocate_workspace();
}

especia::D_Decompose::~D_Decompose() {
}

void especia::D_Decompose::operator()(const R_type A[], R_type Z[],
                                      R_type w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[m * m], Z);

    Z_type info = 0;
    LAPACK_NAME_R_TYPE(syevd)('V', 'U', m, &Z[0], m, w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // To convert from column-major into row-major layout
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::D_Decompose::allocate_workspace() {
    Z_type info = 0;

    LAPACK_NAME_R_TYPE(syevd)('V', 'U', m, 0, m, 0, &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        work.resize(static_cast<N_type>(lwork));
        liwork = iwork[0];
        iwork.resize(static_cast<N_type>(liwork));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::D_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < m; ++i, i0 += m) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += m) {
            swap(A[ij], A[ji]);
        }
    }
}

const string especia::D_Decompose::message_int_err = "especia::D_Decompose() Error: internal error in LAPACK";
const string especia::D_Decompose::message_ill_arg = "especia::D_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::R_Decompose::R_Decompose(N_type n)
        : m(Z_type(n)), work(1), iwork(1), isupp(2 * max(N_type(1), n)), awork(n * n) {
    allocate_workspace();
}

especia::R_Decompose::~R_Decompose() {
}

void especia::R_Decompose::operator()(const R_type A[], R_type Z[],
                                      R_type w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[m * m], &awork[0]);

    Z_type info = 0;
    Z_type mvec = 0;

    LAPACK_NAME_R_TYPE(syevr)('V', 'A', 'U', m, &awork[0], m, 0.0, 0.0, 0, 0, safe, mvec, w, Z, m,
                              &isupp[0], &work[0], lwork, &iwork[0], liwork,
                              info);

    if (info == 0) {
        // To convert from column-major into row-major layout
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::R_Decompose::allocate_workspace() {
    Z_type info = 0;
    Z_type mvec = 0;

    LAPACK_NAME_R_TYPE(syevr)('V', 'A', 'U', m, 0, m, 0.0, 0.0, 0, 0, 0.0, mvec, 0, 0, m,
                              &isupp[0], &work[0], -1, &iwork[0], -1,
                              info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        work.resize(static_cast<N_type>(lwork));
        liwork = iwork[0];
        iwork.resize(static_cast<N_type>(liwork));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::R_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < m; ++i, i0 += m) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += m) {
            swap(A[ij], A[ji]);
        }
    }
}

const string especia::R_Decompose::message_int_err = "especia::R_Decompose() Error: internal error in LAPACK";
const string especia::R_Decompose::message_ill_arg = "especia::R_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::X_Decompose::X_Decompose(N_type n)
        : m(Z_type(n)), work(1), iwork(5 * n), ifail(n), awork(n * n) {
    allocate_workspace();
}

especia::X_Decompose::~X_Decompose() {
}

void especia::X_Decompose::operator()(const R_type A[], R_type Z[],
                                      R_type w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[m * m], &awork[0]);

    Z_type info = 0;
    Z_type mvec = 0;

    LAPACK_NAME_R_TYPE(syevx)('V', 'A', 'U', m, &awork[0], m, 0.0, 0.0, 0, 0, 2.0 * safe, mvec, w, Z, m,
                              &work[0], lwork, &iwork[0], &ifail[0],
                              info);

    if (info == 0) {
        // To convert from column-major into row-major layout
        transpose(Z);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::X_Decompose::allocate_workspace() {
    Z_type info = 0;
    Z_type mvec = 0;

    LAPACK_NAME_R_TYPE(syevx)('V', 'A', 'U', m, 0, m, 0.0, 0.0, 0, 0, 0.0, mvec, 0, 0, m,
                              &work[0], -1, &iwork[0], &ifail[0],
                              info);

    if (info == 0) {
        lwork = static_cast<Z_type>(work[0]);
        work.resize(static_cast<N_type>(lwork));
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::X_Decompose::transpose(R_type A[]) const {
    for (Z_type i = 0, i0 = 0; i < m; ++i, i0 += m) {
        for (Z_type j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += m) {
            swap(A[ij], A[ji]);
        }
    }
}

const string especia::X_Decompose::message_int_err = "especia::X_Decompose() Error: internal error in LAPACK";
const string especia::X_Decompose::message_ill_arg = "especia::X_Decompose() Error: illegal argument(s) in call to LAPACK";

#undef LAPACK_NAME_R_TYPE
#undef LAPACK_NAME_SINGLE
#undef LAPACK_NAME_DOUBLE
