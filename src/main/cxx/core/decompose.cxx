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
using std::valarray;

using especia::Real;
using especia::Integer;


#define LAPACK_NAME_DOUBLE(x) d##x##_
#define LAPACK_NAME_SINGLE(x) s##x##_
#define LAPACK_NAME_R_TYPE(x) LAPACK_NAME_DOUBLE(x)

extern "C" {
/**
 * Interface to LAPACK routine @c [DS]LAMCH.
 */
Real LAPACK_NAME_R_TYPE(lamch)(const char &cmach);

/**
 * Interface to LAPACK routine @c [DS]SYEVD.
 */
void LAPACK_NAME_R_TYPE(syevd)(const char &job,
                               const char &uplo,
                               const Integer &n,
                               Real A[],
                               const Integer &lda,
                               Real w[],
                               Real work[],
                               const Integer &lwork,
                               Integer iwork[],
                               const Integer &liwork,
                               Integer &info);

/**
 * Interface to LAPACK routine @c [DS]SYEVR.
 */
void LAPACK_NAME_R_TYPE(syevr)(const char &job,
                               const char &range,
                               const char &uplo,
                               const Integer &n,
                               Real A[],
                               const Integer &lda,
                               const Real &vl,
                               const Real &vu,
                               const Integer &il,
                               const Integer &iu,
                               const Real &abstol,
                               Integer &m,
                               Real w[],
                               Real Z[],
                               const Integer &ldz,
                               Integer isupp[],
                               Real work[],
                               const Integer &lwork,
                               Integer iwork[],
                               const Integer &liwork,
                               Integer &info);

/**
 * Interface to LAPACK routine @c [DS]SYEVX.
 */
void LAPACK_NAME_R_TYPE(syevx)(const char &job,
                               const char &range,
                               const char &uplo,
                               const Integer &n,
                               Real A[],
                               const Integer &lda,
                               const Real &vl,
                               const Real &vu,
                               const Integer &il,
                               const Integer &iu,
                               const Real &abstol,
                               Integer &m,
                               Real w[],
                               Real Z[],
                               const Integer &ldz,
                               Real work[],
                               const Integer &lwork,
                               Integer iwork[],
                               Integer ifail[],
                               Integer &info);
}

/**
 * The LAPACK job parameter (here: compute eigenvalues and eigenvectors).
 */
static const char job = 'V';

/**
 * The LAPACK range parameter (here: compute all eigenvalues and eigenvectors).
 */
static const char range = 'A';

/**
 * The LAPACK matrix store parameter (here: use the upper triangular part)
 */
static const char uplo = 'U';

/**
 * The LAPACK lower range limit (here: not used).
 */
static const Real vl = 0.0;

/**
 * The LAPACK upper range limit (here: not used).
 */
static const Real vu = 0.0;

/**
 * The LAPACK lower range index (here: not used).
 */
static const Integer il = 0;

/**
 * The LAPACK upper range index (here: not used).
 */
static const Integer iu = 0;


especia::D_Decompose::D_Decompose(Natural m)
        : n(Integer(m)), work(), iwork() {
    lapack_inquire(n, lwork, liwork);

    work.resize((size_t) lwork);
    iwork.resize((size_t) liwork);
}

especia::D_Decompose::~D_Decompose() {
}

void especia::D_Decompose::operator()(const Real A[],
                                      Real Z[],
                                      Real w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[n * n], Z);

    lapack_do(Z, w);
    transpose(n, Z);
}

void especia::D_Decompose::lapack_do(Real Z[], Real w[]) const {
    Integer info = 0;

    LAPACK_NAME_R_TYPE(syevd)(job, uplo, n, &Z[0], n, w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        // ok
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

void especia::D_Decompose::lapack_inquire(Integer n, Integer &lwork, Integer &liwork) {
    Integer info;
    Real work;

    LAPACK_NAME_R_TYPE(syevd)(job, uplo, n, 0, n, 0, &work, -1, &liwork, -1, info);

    if (info == 0) {
        lwork = static_cast<Integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const string especia::D_Decompose::message_int_err = "especia::D_Decompose() Error: internal error in LAPACK";
const string especia::D_Decompose::message_ill_arg = "especia::D_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::R_Decompose::R_Decompose(Natural m)
        : n(Integer(m)), work(), iwork(), isupp(2 * max<Natural>(1, m)), awork(m * m) {
    lapack_inquire(n, lwork, liwork);

    work.resize((size_t) lwork);
    iwork.resize((size_t) liwork);
}

especia::R_Decompose::~R_Decompose() {
}

void especia::R_Decompose::operator()(const Real A[],
                                      Real Z[],
                                      Real w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[n * n], &awork[0]);

    lapack_do(Z, w);
    transpose(n, Z);
}

void especia::R_Decompose::lapack_do(Real Z[], Real w[]) const {
    Integer m = 0;
    Integer info = 0;

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

void especia::R_Decompose::lapack_inquire(Integer n, Integer &lwork, Integer &liwork) {
    Integer info;
    Integer m;
    Real work;

    LAPACK_NAME_R_TYPE(syevr)(job, range, uplo, n, 0, n, vl, vu, il, iu, abstol, m, 0, 0, n,
                              0, &work, -1, &liwork, -1,
                              info);

    if (info == 0) {
        lwork = static_cast<Integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const Real especia::R_Decompose::abstol = LAPACK_NAME_R_TYPE(lamch)('S');
const string especia::R_Decompose::message_int_err = "especia::R_Decompose() Error: internal error in LAPACK";
const string especia::R_Decompose::message_ill_arg = "especia::R_Decompose() Error: illegal argument(s) in call to LAPACK";


especia::X_Decompose::X_Decompose(Natural m)
        : n(Integer(m)), work(), iwork(5 * m), ifail(m), awork(m * m) {
    lapack_inquire(n, lwork);

    work.resize((size_t) lwork);
}

especia::X_Decompose::~X_Decompose() {
}

void especia::X_Decompose::operator()(const Real A[],
                                      Real Z[],
                                      Real w[]) const throw(invalid_argument, runtime_error) {
    copy(&A[0], &A[n * n], &awork[0]);

    lapack_do(Z, w);
    transpose(n, Z);
}

void especia::X_Decompose::lapack_do(Real Z[], Real w[]) const {
    Integer m = 0;
    Integer info = 0;

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

void especia::X_Decompose::lapack_inquire(Integer n, Integer &lwork) {
    Integer info;
    Integer m;
    Real work;

    LAPACK_NAME_R_TYPE(syevx)(job, range, uplo, n, 0, n, vl, vu, il, iu, abstol, m, 0, 0, n,
                              &work, -1, 0, 0,
                              info);

    if (info == 0) {
        lwork = static_cast<Integer>(work);
    } else if (info > 0) {
        throw runtime_error(message_int_err);
    } else {
        throw invalid_argument(message_ill_arg);
    }
}

const Real especia::X_Decompose::abstol = Real(2) * LAPACK_NAME_R_TYPE(lamch)('S');
const string especia::X_Decompose::message_int_err = "especia::X_Decompose() Error: internal error in LAPACK";
const string especia::X_Decompose::message_ill_arg = "especia::X_Decompose() Error: illegal argument(s) in call to LAPACK";

#undef LAPACK_NAME_R_TYPE
#undef LAPACK_NAME_SINGLE
#undef LAPACK_NAME_DOUBLE
