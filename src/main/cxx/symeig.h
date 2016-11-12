// Symmetric eigenproblem solvers calling LAPACK
// Copyright (c) 2016, Ralf Quast
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RQ_SYMEIG_H
#define RQ_SYMEIG_H

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
        size_t n) throw (std::runtime_error);

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
        size_t n) throw (std::runtime_error);

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
        size_t n) throw (std::runtime_error);

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
