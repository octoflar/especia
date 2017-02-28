// Symmetric eigenproblem solvers calling LAPACK
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
#ifndef ESPECIA_SYMEIG_H
#define ESPECIA_SYMEIG_H

#include <cstddef>
#include <stdexcept>
#include <valarray>

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
         * @param n The problem dimension.
         */
        D_Decompose(size_t n = 0);

        /**
         * Destructor.
         */
        ~D_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param n[in] The problem dimension.
         * @param A[in] The symmetric matrix (row-major, lower triangular).
         * @param Z[out] The transformation matrix (row-major).
         * @param w[out] The eigenvalues.
         */
        void operator()(size_t n, const double A[], double Z[], double w[]) throw(std::runtime_error);

    private:
        void resize_workspace(size_t n = 0);

        void transpose(double A[]) const;

        const char job;
        const char uplo;

        int n;
        int info;

        int lwork;
        int liwork;
        std::valarray<double> work;
        std::valarray<int> iwork;

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
         * @param n The problem dimension.
         */
        R_Decompose(size_t n = 0);

        /**
         * Destructor.
         */
        ~R_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param n[in] The problem dimension.
         * @param A[in] The symmetric matrix (row-major, lower triangular).
         * @param Z[out] The transformation matrix (row-major).
         * @param w[out] The eigenvalues.
         */
        void operator()(size_t n, const double A[], double Z[], double w[]) throw(std::runtime_error);

    private:
        void resize_workspace(size_t n = 0);

        void transpose(double A[]) const;

        const char job;
        const char range;
        const char uplo;

        int m;
        int n;
        int info;

        std::valarray<int> isupp;

        int lwork;
        int liwork;
        std::valarray<double> work;
        std::valarray<int> iwork;

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
         * @param n The problem dimension.
         */
        X_Decompose(size_t n = 0);

        /**
         * Destructor.
         */
        ~X_Decompose();

        /**
         * Solves a symmetric eigenproblem.
         *
         * @param n[in] The problem dimension.
         * @param A[in] The symmetric matrix (row-major, lower triangular).
         * @param Z[out] The transformation matrix (row-major).
         * @param w[out] The eigenvalues.
         */
        void operator()(size_t n, const double A[], double Z[], double w[]) throw(std::runtime_error);

    private:
        void resize_workspace(size_t n = 0);

        void transpose(double A[]) const;

        const char job;
        const char range;
        const char uplo;

        int m;
        int n;
        int info;

        int lwork;
        std::valarray<double> work;
        std::valarray<int> iwork;
        std::valarray<int> ifail;

        static const char int_err[];
        static const char ill_arg[];
    };

    /**
     * The selected algorithm to solve symmetric eigenproblems.
     */
    typedef R_Decompose Decompose;

}

#endif // ESPECIA_SYMEIG_H