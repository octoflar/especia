// Symmetric eigenproblem solver calling the LAPACK routine DSYEVR
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
#include <algorithm>
#include "decompose.h"

#if !defined(F77NAME)
#define F77NAME(x) x##_
#endif

using std::copy;
using std::max;
using std::runtime_error;
using std::swap;
using std::valarray;

// Interface to LAPACK eigenvalue routines (version 3.0)
extern "C" {
double F77NAME(dlamch)(const char &cmach);

void F77NAME(dsyevd)(const char &job,
                     const char &uplo,
                     const int &n,
                     double A[], const int &lda,
                     double w[],
                     double work[], const int &lwork,
                     int iwork[], const int &liwork,
                     int &info);

void F77NAME(dsyevr)(const char &job,
                     const char &range,
                     const char &uplo,
                     const int &n,
                     double A[], const int &lda,
                     const double &vl, const double &vu,
                     const int &il, const int &iu,
                     const double &abstol,
                     int &m,
                     double w[],
                     double Z[], const int &ldz,
                     int isupp[],
                     double work[], const int &lwork,
                     int iwork[], const int &liwork,
                     int &info);

void F77NAME(dsyevx)(const char &job,
                     const char &range,
                     const char &uplo,
                     const int &n,
                     double A[], const int &lda,
                     const double &vl, const double &vu,
                     const int &il, const int &iu,
                     const double &abstol,
                     int &m,
                     double w[],
                     double Z[], const int &ldz,
                     double work[], const int &lwork,
                     int iwork[],
                     int ifail[], int &info);
}

const double safe_minimum = F77NAME(dlamch)('s');

const char especia::D_Decompose::int_err[] = "especia::D_Decompose(): Error: internal error in LAPACK routine DSYEVD";
const char especia::D_Decompose::ill_arg[] = "especia::D_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVD";

especia::D_Decompose::D_Decompose(size_t n)
        : job('V'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::D_Decompose::~D_Decompose() {
}

void especia::D_Decompose::operator()(const double A[], double Z[], double w[], size_t k) throw(runtime_error) {
    copy(&A[0], &A[k * k], Z);

    if (k != n)
        resize_workspace(k);

    // regular call
    F77NAME(dsyevd)(job, uplo, n, &Z[0], max(1, n), w, &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        transpose(Z);
            // transform from column-major into row-major layout
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::D_Decompose::resize_workspace(size_t k) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
    n = k;
#pragma clang diagnostic pop

    // workspace query
    F77NAME(dsyevd)(job, uplo, n, 0, max(1, n), 0, &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<int>(work[0]);
        liwork = iwork[0];

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        work.resize(lwork);
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        iwork.resize(liwork);
#pragma clang diagnostic pop
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::D_Decompose::transpose(double A[]) const {
    for (int i = 0, i0 = 0; i < n; ++i, i0 += n)
        for (int j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n)
            swap(A[ij], A[ji]);
}


const char especia::R_Decompose::int_err[] = "especia::R_Decompose(): Error: internal error in LAPACK routine DSYEVR";
const char especia::R_Decompose::ill_arg[] = "especia::R_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVR";

especia::R_Decompose::R_Decompose(size_t n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(1) {
    resize_workspace(n);
}

especia::R_Decompose::~R_Decompose() {
}

void especia::R_Decompose::operator()(const double A[], double Z[], double w[], size_t k) throw(runtime_error) {
    valarray<double> C(A, k * k);

    if (k != n)
        resize_workspace(k);

    // regular call
    F77NAME(dsyevr)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, safe_minimum,
                    m, w, Z, max(1, n), &isupp[0], &work[0], lwork, &iwork[0], liwork, info);

    if (info == 0) {
        transpose(Z);
            // transform from column-major into row-major layout
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::R_Decompose::resize_workspace(size_t k) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
    n = k;
#pragma clang diagnostic pop

    // workspace query
    F77NAME(dsyevr)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, safe_minimum,
                    m, 0, 0, max(1, n), &isupp[0], &work[0], -1, &iwork[0], -1, info);

    if (info == 0) {
        lwork = static_cast<int>(work[0]);
        liwork = iwork[0];

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        work.resize(lwork);
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        iwork.resize(liwork);
#pragma clang diagnostic pop

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        isupp.resize(2 * max(1, n));
#pragma clang diagnostic pop
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::R_Decompose::transpose(double A[]) const {
    for (int i = 0, i0 = 0; i < n; ++i, i0 += n)
        for (int j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n)
            swap(A[ij], A[ji]);
}


const char especia::X_Decompose::int_err[] = "especia::X_Decompose(): Error: internal error in LAPACK routine DSYEVX";
const char especia::X_Decompose::ill_arg[] = "especia::X_Decompose(): Error: illegal argument(s) in call to LAPACK routine DSYEVX";

especia::X_Decompose::X_Decompose(size_t n)
        : job('V'), range('A'), uplo('U'), work(1), iwork(), ifail() {
    resize_workspace(n);
}

especia::X_Decompose::~X_Decompose() {
}

void especia::X_Decompose::operator()(const double A[], double Z[], double w[], size_t k) throw(runtime_error) {
    valarray<double> C(A, k * k);

    if (k != n)
        resize_workspace(k);

    // regular call
    F77NAME(dsyevx)(job, range, uplo, n, &C[0], max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum,
                    m, w, Z, max(1, n), &work[0], lwork, &iwork[0], &ifail[0], info);

    if (info == 0) {
        transpose(Z);
            // transform from column-major into row-major layout
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::X_Decompose::resize_workspace(size_t k) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
    n = k;
#pragma clang diagnostic pop

    // workspace query
    F77NAME(dsyevx)(job, range, uplo, n, 0, max(1, n), 0.0, 0.0, 0, 0, 2.0 * safe_minimum,
                    m, 0, 0, max(1, n), &work[0], -1, &iwork[0], &ifail[0], info);

    if (info == 0) {
        lwork = static_cast<int>(work[0]);

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        work.resize(max(1, lwork));
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        iwork.resize(5 * n);
#pragma clang diagnostic pop
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
        ifail.resize(n);
#pragma clang diagnostic pop
    } else if (info > 0)
        throw runtime_error(int_err);
    else
        throw runtime_error(ill_arg);
}

void especia::X_Decompose::transpose(double A[]) const {
    for (int i = 0, i0 = 0; i < n; ++i, i0 += n)
        for (int j = 0, ij = i0, ji = i; j < i; ++j, ++ij, ji += n)
            swap(A[ij], A[ji]);
}
