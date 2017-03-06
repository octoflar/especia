/// @file decompose.h
/// Symmetric eigenproblem solvers calling LAPACK routines.
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
#ifndef ESPECIA_SYMEIG_H
#define ESPECIA_SYMEIG_H

#include <cstddef>
#include <stdexcept>
#include <valarray>

#include "base.h"

namespace especia {

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * DSYEVD (divide and conquer).
     */
    class D_Decompose {
    public:
        /**
         * Constructs a new instance of this class.
         *
         * @param[in] n The problem dimension.
         */
        D_Decompose(Nnum_t n = 0);

        /**
         * Destructor.
         */
        ~D_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] n The problem dimension.
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues.
         */
        void operator()(Nnum_t n, const Real_t A[], Real_t Z[], Real_t w[]) throw(std::runtime_error);

    private:
        void resize_workspace(Nnum_t n = 0);

        void transpose(Real_t A[]) const;

        const char job;
        const char uplo;

        Znum_t n;
        Znum_t info;

        Znum_t lwork;
        Znum_t liwork;
        std::valarray<Real_t> work;
        std::valarray<Znum_t> iwork;

        static const char int_err[];
        static const char ill_arg[];
    };

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * DSYEVR (relatively robust representations).
     */
    class R_Decompose {
    public:
        /**
         * Constructs a new instance of this class.
         *
         * @param[in] n The problem dimension.
         */
        R_Decompose(Nnum_t n = 0);

        /**
         * Destructor.
         */
        ~R_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] n The problem dimension.
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues.
         */
        void operator()(Nnum_t n, const Real_t A[], Real_t Z[], Real_t w[]) throw(std::runtime_error);

    private:
        void resize_workspace(Nnum_t n = 0);

        void transpose(Real_t A[]) const;

        const char job;
        const char range;
        const char uplo;

        Znum_t m;
        Znum_t n;
        Znum_t info;

        std::valarray<Znum_t> isupp;

        Znum_t lwork;
        Znum_t liwork;
        std::valarray<Real_t> work;
        std::valarray<Znum_t> iwork;

        static const char int_err[];
        static const char ill_arg[];
    };

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * DSYEVX (inverse iteration).
     */
    class X_Decompose {
    public:
        /**
         * Constructs a new instance of this class.
         *
         * @param[in] n The problem dimension.
         */
        X_Decompose(Nnum_t n = 0);

        /**
         * Destructor.
         */
        ~X_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] n The problem dimension.
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues.
         */
        void operator()(Nnum_t n, const Real_t A[], Real_t Z[], Real_t w[]) throw(std::runtime_error);

    private:
        void resize_workspace(Nnum_t n = 0);

        void transpose(Real_t A[]) const;

        const char job;
        const char range;
        const char uplo;

        Znum_t m;
        Znum_t n;
        Znum_t info;

        Znum_t lwork;
        std::valarray<Real_t> work;
        std::valarray<Znum_t> iwork;
        std::valarray<Znum_t> ifail;

        static const char int_err[];
        static const char ill_arg[];
    };

    /**
     * The selected algorithm to solve symmetric eigenproblems.
     */
    typedef R_Decompose Decompose;

}

#endif // ESPECIA_SYMEIG_H
