/// @file random_test.cxx
/// Unit tests
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include "../../../main/cxx/core/random.h"
#include "../unittest.h"

using especia::Melg19937_64;
using especia::Mt19937_32;
using especia::Mt19937_64;
using especia::Pcg_32;

class Rng_Test : public Unit_Test
{
private:
  void
  test_melg19937_64 ()
  {
    using especia::natural;
    using especia::word64;

    const word64 seeds[] = { 0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL };
    const Melg19937_64 melg (4, seeds);

    // <https://github.com/sharase/melg-64>
    assert_equals (16675511042081433281ULL, melg.rand (),
                   "test MELG-19937-64 (0)");
    assert_equals (8489326016911908102ULL, melg.rand (),
                   "test MELG-19937-64 (1)");
    assert_equals (16071362722047509693ULL, melg.rand (),
                   "test MELG-19937-64 (2)");
    assert_equals (11631833934008589069ULL, melg.rand (),
                   "test MELG-19937-64 (3)");
    assert_equals (3308423691540511443ULL, melg.rand (),
                   "test MELG-19937-64 (4)");
    assert_equals (12463994900921303743ULL, melg.rand (),
                   "test MELG-19937-64 (5)");

    for (especia::natural i = 6; i < 999; i++)
      {
        melg.rand ();
      }

    assert_equals (13711744326396256691ULL, melg.rand (),
                   "test MELG-19937-64 (6)");
  }

  void
  test_mt19937_32 ()
  {
    using especia::word64;

    const word64 seeds[] = { 0x123ULL, 0x234ULL, 0x345ULL, 0x456ULL };
    const Mt19937_32 mt (4, seeds);

    assert_equals (1067595299ULL, mt.rand (), "test MT-19937-32 (0)");
    assert_equals (955945823ULL, mt.rand (), "test MT-19937-32 (1)");
    assert_equals (477289528ULL, mt.rand (), "test MT-19937-32 (2)");
    assert_equals (4107218783ULL, mt.rand (), "test MT-19937-32 (3)");
    assert_equals (4228976476ULL, mt.rand (), "test MT-19937-32 (4)");
    assert_equals (3344332714ULL, mt.rand (), "test MT-19937-32 (5)");
  }

  void
  test_mt19937_64 ()
  {
    using especia::word64;

    const word64 seeds[] = { 0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL };
    const Mt19937_64 mt (4, seeds);

    assert_equals (7266447313870364031ULL, mt.rand (), "test MT-19937-64 (0)");
    assert_equals (4946485549665804864ULL, mt.rand (), "test MT-19937-64 (1)");
    assert_equals (16945909448695747420ULL, mt.rand (),
                   "test MT-19937-64 (2)");
    assert_equals (16394063075524226720ULL, mt.rand (),
                   "test MT-19937-64 (3)");
    assert_equals (4873882236456199058ULL, mt.rand (), "test MT-19937-64 (4)");
    assert_equals (14877448043947020171ULL, mt.rand (),
                   "test MT-19937-64 (5)");
  }

  void
  test_pcg ()
  {
    const Pcg_32 pcg (42ULL, 54ULL);

    assert_equals (0xa15c02b7UL, pcg.rand (), "test PCG-XSH-RR-64-32 (0)");
    assert_equals (0x7b47f409UL, pcg.rand (), "test PCG-XSH-RR-64-32 (1)");
    assert_equals (0xba1d3330UL, pcg.rand (), "test PCG-XSH-RR-64-32 (2)");
    assert_equals (0x83d2f293UL, pcg.rand (), "test PCG-XSH-RR-64-32 (3)");
    assert_equals (0xbfa4784bUL, pcg.rand (), "test PCG-XSH-RR-64-32 (4)");
    assert_equals (0xcbed606eUL, pcg.rand (), "test PCG-XSH-RR-64-32 (5)");
  }

  void
  run_all () override
  {
    run (this, &Rng_Test::test_melg19937_64);
    run (this, &Rng_Test::test_mt19937_32);
    run (this, &Rng_Test::test_mt19937_64);
    run (this, &Rng_Test::test_pcg);
  }
};

int
main ()
{
  return Rng_Test ().run_testsuite ();
}
