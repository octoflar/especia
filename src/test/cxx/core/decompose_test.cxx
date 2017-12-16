/// @file decompose_test.cxx
/// Unit tests
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
#include "../../../main/cxx/core/decompose.h"
#include "../unittest.h"

using especia::Decompose;
using especia::natural;
using especia::real;

class Decompose_Test : public Unit_Test {
private:

    void test_decompose_diagonal_matrix() {
        const natural n = 3;
        const Decompose decompose(n);

        const real A[n * n] = { real(1), real(0), real(0),
                                real(0), real(2), real(0),
                                real(0), real(0), real(3) };

        real Z[n * n];
        real w[n];

        decompose(A, Z, w);

        assert_equals(real(1), Z[0], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[1], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[2], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[3], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(1), Z[4], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[5], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[6], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(0), Z[7], real(0), "decompose diagonal matrix (Z)");
        assert_equals(real(1), Z[8], real(0), "decompose diagonal matrix (Z)");

        assert_equals(real(1), w[0], real(0), "decompose diagonal matrix (w)");
        assert_equals(real(2), w[1], real(0), "decompose diagonal matrix (w)");
        assert_equals(real(3), w[2], real(0), "decompose diagonal matrix (w)");
    }

    void test_decompose_symmetric_matrix() {
        const natural n = 3;
        const Decompose decompose(n);

        const real A[n * n] = { real(1), real(2), real(3),
                                real(2), real(4), real(5),
                                real(3), real(5), real(6) };

        real Z[n * n];
        real w[n];

        decompose(A, Z, w);

        assert_equals(real(-0.736976), Z[0], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real(-0.327985), Z[1], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real( 0.591009), Z[2], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real( 0.591009), Z[3], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real(-0.736976), Z[4], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real( 0.327985), Z[5], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real(-0.327985), Z[6], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real(-0.591009), Z[7], real(1.0E-06), "decompose symmetric matrix (Z)");
        assert_equals(real(-0.736976), Z[8], real(1.0E-06), "decompose symmetric matrix (Z)");

        assert_equals(real(-0.515729), w[0], real(1.0E-06), "decompose symmetric matrix (w)");
        assert_equals(real( 0.170915), w[1], real(1.0E-06), "decompose symmetric matrix (w)");
        assert_equals(real( 11.34480), w[2], real(1.0E-04), "decompose symmetric matrix (w)");
    }

    void run_all() override {
        run(this, &Decompose_Test::test_decompose_diagonal_matrix);
        run(this, &Decompose_Test::test_decompose_symmetric_matrix);
    }
};


int main() {
    return Decompose_Test().run_testsuite();
}
