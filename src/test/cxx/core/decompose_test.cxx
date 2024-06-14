/// @file decompose_test.cxx
/// Unit tests
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include "../../../main/cxx/core/decompose.h"
#include "../unittest.h"

using especia::natural;
using especia::real;

class Decompose_Test : public Unit_Test
{
private:
  void
  test_decompose_diagonal_matrix_D ()
  {
    using especia::D_Decompose;

    const natural n = 3;
    const D_Decompose decompose (n);

    const real A[n * n] = { real (1), real (0), real (0), real (0), real (2),
                            real (0), real (0), real (0), real (3) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    assert_equals (real (1), Z[0], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[1], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[2], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[3], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[4], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[5], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[6], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[7], real (0),
                   "D decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[8], real (0),
                   "D decompose diagonal matrix (Z)");

    assert_equals (real (1), w[0], real (0),
                   "D decompose diagonal matrix (w)");
    assert_equals (real (2), w[1], real (0),
                   "D decompose diagonal matrix (w)");
    assert_equals (real (3), w[2], real (0),
                   "D decompose diagonal matrix (w)");
  }

  void
  test_decompose_diagonal_matrix_R ()
  {
    using especia::R_Decompose;

    const natural n = 3;
    const R_Decompose decompose (n);

    const real A[n * n] = { real (1), real (0), real (0), real (0), real (2),
                            real (0), real (0), real (0), real (3) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    assert_equals (real (1), Z[0], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[1], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[2], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[3], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[4], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[5], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[6], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[7], real (0),
                   "R decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[8], real (0),
                   "R decompose diagonal matrix (Z)");

    assert_equals (real (1), w[0], real (0),
                   "R decompose diagonal matrix (w)");
    assert_equals (real (2), w[1], real (0),
                   "R decompose diagonal matrix (w)");
    assert_equals (real (3), w[2], real (0),
                   "R decompose diagonal matrix (w)");
  }

  void
  test_decompose_diagonal_matrix_X ()
  {
    using especia::X_Decompose;

    const natural n = 3;
    const X_Decompose decompose (n);

    const real A[n * n] = { real (1), real (0), real (0), real (0), real (2),
                            real (0), real (0), real (0), real (3) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    assert_equals (real (1), Z[0], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[1], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[2], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[3], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[4], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[5], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[6], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (0), Z[7], real (0),
                   "X decompose diagonal matrix (Z)");
    assert_equals (real (1), Z[8], real (0),
                   "X decompose diagonal matrix (Z)");

    assert_equals (real (1), w[0], real (0),
                   "X decompose diagonal matrix (w)");
    assert_equals (real (2), w[1], real (0),
                   "X decompose diagonal matrix (w)");
    assert_equals (real (3), w[2], real (0),
                   "X decompose diagonal matrix (w)");
  }

  void
  test_decompose_symmetric_matrix_D ()
  {
    using especia::D_Decompose;

    const natural n = 3;
    const D_Decompose decompose (n);

    const real A[n * n] = { real (1), real (2), real (3), real (2), real (4),
                            real (5), real (3), real (5), real (6) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    // the eigenvectors have different signs than those provided by Wolfram
    // Alpha
    assert_equals (real (0.736976), Z[0], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[1], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (-0.591009), Z[2], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[3], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (-0.736976), Z[4], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[5], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[6], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[7], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");
    assert_equals (real (0.736976), Z[8], real (1.0E-06),
                   "D decompose symmetric matrix (Z)");

    assert_equals (real (-0.515729), w[0], real (1.0E-06),
                   "D decompose symmetric matrix (w)");
    assert_equals (real (0.170915), w[1], real (1.0E-06),
                   "D decompose symmetric matrix (w)");
    assert_equals (real (11.34480), w[2], real (1.0E-04),
                   "D decompose symmetric matrix (w)");
  }

  void
  test_decompose_symmetric_matrix_R ()
  {
    using especia::R_Decompose;

    const natural n = 3;
    const R_Decompose decompose (n);

    const real A[n * n] = { real (1), real (2), real (3), real (2), real (4),
                            real (5), real (3), real (5), real (6) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    // <https://www.wolframalpha.com/input/?i=eigensystem%5BN%5B%7B%7B1,2,3%7D,%7B2,4,5%7D,%7B3,5,6%7D%7D%5D%5D>
    assert_equals (real (-0.736976), Z[0], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (-0.327985), Z[1], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[2], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[3], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (-0.736976), Z[4], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[5], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (-0.327985), Z[6], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (-0.591009), Z[7], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");
    assert_equals (real (-0.736976), Z[8], real (1.0E-06),
                   "R decompose symmetric matrix (Z)");

    assert_equals (real (-0.515729), w[0], real (1.0E-06),
                   "R decompose symmetric matrix (w)");
    assert_equals (real (0.170915), w[1], real (1.0E-06),
                   "R decompose symmetric matrix (w)");
    assert_equals (real (11.34480), w[2], real (1.0E-04),
                   "R decompose symmetric matrix (w)");
  }

  void
  test_decompose_symmetric_matrix_X ()
  {
    using especia::X_Decompose;

    const natural n = 3;
    const X_Decompose decompose (n);

    const real A[n * n] = { real (1), real (2), real (3), real (2), real (4),
                            real (5), real (3), real (5), real (6) };

    real Z[n * n];
    real w[n];

    decompose (A, Z, w);

    // the eigenvectors have different signs than those provided by Wolfram
    // Alpha
    assert_equals (real (-0.736976), Z[0], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (-0.327985), Z[1], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[2], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[3], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (-0.736976), Z[4], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[5], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.327985), Z[6], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.591009), Z[7], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");
    assert_equals (real (0.736976), Z[8], real (1.0E-06),
                   "X decompose symmetric matrix (Z)");

    assert_equals (real (-0.515729), w[0], real (1.0E-06),
                   "X decompose symmetric matrix (w)");
    assert_equals (real (0.170915), w[1], real (1.0E-06),
                   "X decompose symmetric matrix (w)");
    assert_equals (real (11.34480), w[2], real (1.0E-04),
                   "X decompose symmetric matrix (w)");
  }

  void
  run_all () override
  {
    run (this, &Decompose_Test::test_decompose_diagonal_matrix_D);
    run (this, &Decompose_Test::test_decompose_diagonal_matrix_R);
    run (this, &Decompose_Test::test_decompose_diagonal_matrix_X);
    run (this, &Decompose_Test::test_decompose_symmetric_matrix_D);
    run (this, &Decompose_Test::test_decompose_symmetric_matrix_R);
    run (this, &Decompose_Test::test_decompose_symmetric_matrix_X);
  }
};

int
main ()
{
  return Decompose_Test ().run_testsuite ();
}
