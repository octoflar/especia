/// @file random_test.cxx
/// Unit tests
/// Copyright (c) 2020 Ralf Quast
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
#include "../../../main/cxx/core/random.h"
#include "../unittest.h"

using especia::Mt19937_32;
using especia::Mt19937_64;
using especia::Pcg_32;
using especia::word64;

class Rng_Test : public Unit_Test {
private:

    void test_mt19937_32() {
      const word64 seeds[] = {0x123ull, 0x234ull, 0x345ull, 0x456ull};
      const Mt19937_32 mt(4, seeds);

      assert_equals(1067595299ull, mt.rand(), "test MT-19937-32 (0)");
      assert_equals(955945823ull, mt.rand(), "test MT-19937-32 (1)");
      assert_equals(477289528ull, mt.rand(), "test MT-19937-32 (2)");
      assert_equals(4107218783ull, mt.rand(), "test MT-19937-32 (3)");
      assert_equals(4228976476ull, mt.rand(), "test MT-19937-32 (4)");
      assert_equals(3344332714ull, mt.rand(), "test MT-19937-32 (5)");
    }
  
    void test_mt19937_64() {
      const word64 seeds[] = {0x12345ull, 0x23456ull, 0x34567ull, 0x45678ull};
      const Mt19937_64 mt(4, seeds);

      assert_equals(7266447313870364031ull, mt.rand(), "test MT-19937-64 (0)");
      assert_equals(4946485549665804864ull, mt.rand(), "test MT-19937-64 (1)");
      assert_equals(16945909448695747420ull, mt.rand(), "test MT-19937-64 (2)");
      assert_equals(16394063075524226720ull, mt.rand(), "test MT-19937-64 (3)");
      assert_equals(4873882236456199058ull, mt.rand(), "test MT-19937-64 (4)");
      assert_equals(14877448043947020171ull, mt.rand(), "test MT-19937-64 (5)");
    }

    void test_pcg() {
      const Pcg_32 pcg(42ull, 54ull);
      
      assert_equals(0xa15c02b7ul, pcg.rand(), "test PCG-XSH-RR-64-32 (0)");
      assert_equals(0x7b47f409ul, pcg.rand(), "test PCG-XSH-RR-64-32 (1)");
      assert_equals(0xba1d3330ul, pcg.rand(), "test PCG-XSH-RR-64-32 (2)");
      assert_equals(0x83d2f293ul, pcg.rand(), "test PCG-XSH-RR-64-32 (3)");
      assert_equals(0xbfa4784bul, pcg.rand(), "test PCG-XSH-RR-64-32 (4)");
      assert_equals(0xcbed606eul, pcg.rand(), "test PCG-XSH-RR-64-32 (5)");
    }
  
    
    void run_all() override {
        run(this, &Rng_Test::test_mt19937_32);
        run(this, &Rng_Test::test_mt19937_64);
        run(this, &Rng_Test::test_pcg);
    }
};


int main() {
    return Rng_Test().run_testsuite();
}
