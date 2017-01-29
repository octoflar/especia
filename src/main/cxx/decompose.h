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
    // Function-like class for solving symmetric eigenproblems
    class D_Decompose; // divide and conquer
    class R_Decompose; // relatively robust representations
    class X_Decompose; // inverse iteration
}

// Function-like class for solving symmetric eigenproblems. Calls the
// LAPACK driver routine DSYEVD (divide and conquer).
class especia::D_Decompose {
public:
    D_Decompose(size_t n = 0);

    ~D_Decompose();

    void operator()(const double A[], // symmetric matrix (row-major, lower triangular)
                    double Z[], // transformation matrix (row-major)
                    double w[], // diagonal elements
                    size_t n) throw(std::runtime_error);

private:
    void resize_workspace(size_t n = 0);

    void transpose(double A[]) const;

    char job;
    char uplo;

    int n;
    int info;

    int lwork;
    int liwork;
    std::valarray<double> work;
    std::valarray<int> iwork;

    static const char int_err[];
    static const char ill_arg[];
};

// Function-like class for solving symmetric eigenproblems. Calls the
// LAPACK driver routine DSYEVR (relatively robust representations).
class especia::R_Decompose {
public:
    R_Decompose(size_t n = 0);

    ~R_Decompose();

    void operator()(const double A[], // symmetric matrix (row-major, lower triangular)
                    double Z[], // transformation matrix (row-major)
                    double w[], // diagonal elements
                    size_t n) throw(std::runtime_error);

private:
    void resize_workspace(size_t n = 0);

    void transpose(double A[]) const;

    char job;
    char range;
    char uplo;

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

// Function-like class for solving symmetric eigenproblems. Calls the
// LAPACK driver routine DSYEVX (inverse iteration).
class especia::X_Decompose {
public:
    X_Decompose(size_t n = 0);

    ~X_Decompose();

    void operator()(const double A[], // symmetric matrix (row-major, lower triangular)
                    double Z[], // transformation matrix (row-major)
                    double w[], // diagonal elements
                    size_t n) throw(std::runtime_error);

private:
    void resize_workspace(size_t n = 0);

    void transpose(double A[]) const;

    char job;
    char range;
    char uplo;

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

#endif // ESPECIA_SYMEIG_H
