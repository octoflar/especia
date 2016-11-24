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
#ifndef RQ_SYMEIG_H
#define RQ_SYMEIG_H

#include <cstddef>
#include <stdexcept>
#include <valarray>

namespace RQ {
    // Function-like class for solving symmetric eigenproblems
    class sym_eig_decomp_d; // divide and conquer
    class sym_eig_decomp_r; // relatively robust representations
    class sym_eig_decomp_x; // inverse iteration

    typedef sym_eig_decomp_r sym_eig_decomp;
}

// Function-like class for solving symmetric eigenproblems. Calls the
// LAPACK driver routine DSYEVD (divide and conquer).
class RQ::sym_eig_decomp_d {
public:
    sym_eig_decomp_d(size_t n = 0);

    ~sym_eig_decomp_d();

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
class RQ::sym_eig_decomp_r {
public:
    sym_eig_decomp_r(size_t n = 0);

    ~sym_eig_decomp_r();

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
class RQ::sym_eig_decomp_x {
public:
    sym_eig_decomp_x(size_t n = 0);

    ~sym_eig_decomp_x();

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

#endif // RQ_SYMEIG_H
