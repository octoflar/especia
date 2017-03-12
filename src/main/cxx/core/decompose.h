/// @file decompose.h
/// Symmetric eigenproblem solvers calling LAPACK routines.
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
#ifndef ESPECIA_SYMEIG_H
#define ESPECIA_SYMEIG_H

#include <stdexcept>
#include <string>
#include <valarray>

#include "base.h"

namespace especia {

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * @c [DS]SYEVD (divide and conquer).
     *
     * The divide and conquer algorithm makes very mild assumptions about
     * floating point arithmetic.
     */
    class D_Decompose {
    public:
        /**
         * Constructs a new instance of this class for the problem dimension supplied as argument.
         *
         * @param[in] m The problem dimension.
         */
        D_Decompose(N_type m);

        /**
         * Destructor.
         */
        ~D_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues, in ascending order.
         *
         * @throw invalid_argument when LAPACK was called with illegal arguments.
         * @throw runtime_error when an internal LAPACK error occurred.
         */
        void
        operator()(const R_type A[], R_type Z[], R_type w[]) const throw(std::invalid_argument, std::runtime_error);

    private:
        /**
         * Queries and allocates the optimal workspace.
         */
        void allocate_workspace();

        /**
         * Transposes a matrix. Called to convert from column-major into row-major layout.
         *
         * @param[in,out] A The matrix.
         */
        void transpose(R_type A[]) const;

        /**
         * The problem dimension.
         */
        const Z_type n;

        /**
         * A workspace size.
         */
        Z_type lwork;

        /**
         * A workspace size.
         */
        Z_type liwork;

        /**
         * A workspace array.
         */
        mutable std::valarray<R_type> work;

        /**
         * A workspace array.
         */
        mutable std::valarray<Z_type> iwork;

        static const std::string message_int_err;
        static const std::string message_ill_arg;
    };

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * @c [DS]SYEVR (relatively robust representations).
     *
     * Normal execution may create NaN and infinities and hence may abort
     * due to a floating point exception in environments, which do not
     * handle NaN and infinities in the IEEE standard default manner.
     */
    class R_Decompose {
    public:
        /**
         * Constructs a new instance of this class for the problem dimension supplied as argument.
         *
         * @param[in] m The problem dimension.
         */
        R_Decompose(N_type m);

        /**
         * Destructor.
         */
        ~R_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues, in ascending order.
         *
         * @throw invalid_argument when LAPACK was called with illegal arguments.
         * @throw runtime_error when an internal LAPACK error occurred.
         */
        void
        operator()(const R_type A[], R_type Z[], R_type w[]) const throw(std::invalid_argument, std::runtime_error);

    private:
        /**
         * Queries and allocates the optimal workspace.
         */
        void allocate_workspace();

        /**
         * Transposes a matrix. Called to convert from column-major into row-major layout.
         *
         * @param[in,out] A The matrix.
         */
        void transpose(R_type A[]) const;

        /**
         * The problem dimension.
         */
        const Z_type n;

        /**
         * A workspace size.
         */
        Z_type lwork;

        /**
         * A workspace size.
         */
        Z_type liwork;

        /**
         * A workspace array.
         */
        mutable std::valarray<R_type> work;

        /**
         * A workspace array.
         */
        mutable std::valarray<Z_type> iwork;

        /**
         * A workspace array.
         */
        mutable std::valarray<Z_type> isupp;

        /**
         * A workspace array.
         */
        mutable std::valarray<R_type> awork;

        /**
         * The absolute accuracy of eigenvalues computed. Yields the most accurate results
         * when set to the 'safe minimum'.
         */
        static const R_type abstol;

        static const std::string message_int_err;
        static const std::string message_ill_arg;
    };

    /**
     * Class to solve symmetric eigenproblems. Calls the LAPACK driver routine
     * @c [DS]SYEVX (inverse iteration).
     */
    class X_Decompose {
    public:
        /**
         * Constructs a new instance of this class for the problem dimension supplied as argument.
         *
         * @param[in] m The problem dimension.
         */
        X_Decompose(N_type m = 0);

        /**
         * Destructor.
         */
        ~X_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param[in] A The symmetric matrix (row-major, lower triangular).
         * @param[out] Z The transformation matrix (row-major).
         * @param[out] w The eigenvalues, in ascending order.
         *
         * @throw invalid_argument when LAPACK was called with illegal arguments.
         * @throw runtime_error when an internal LAPACK error occurred.
         */
        void
        operator()(const R_type A[], R_type Z[], R_type w[]) const throw(std::invalid_argument, std::runtime_error);

    private:
        /**
         * Queries and allocates the optimal workspace.
         */
        void allocate_workspace();

        /**
         * Transposes a matrix. Called to convert from column-major into row-major layout.
         *
         * @param[in,out] A The matrix.
         */
        void transpose(R_type A[]) const;

        /**
         * The problem dimension.
         */
        const Z_type n;

        /**
         * A workspace size.
         */
        Z_type lwork;

        /**
         * A workspace array.
         */
        mutable std::valarray<R_type> work;

        /**
         * A workspace array.
         */
        mutable std::valarray<Z_type> iwork;

        /**
         * A workspace array.
         */
        mutable std::valarray<Z_type> ifail;

        /**
         * A workspace array.
         */
        mutable std::valarray<R_type> awork;

        /**
         * The absolute accuracy of eigenvalues computed. Yields the most accurate results
         * when set to twice the 'safe minimum'.
         */
        static const R_type abstol;

        static const std::string message_int_err;
        static const std::string message_ill_arg;
    };

    /**
     * The selected algorithm to solve symmetric eigenproblems.
     */
    typedef R_Decompose Decompose;

}

#endif // ESPECIA_SYMEIG_H
