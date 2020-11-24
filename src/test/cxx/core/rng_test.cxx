/// @file rng_test.cxx
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
#include "../../../main/cxx/core/rng.h"
#include "../unittest.h"

using especia::Mt19937;
using especia::Pcg32;
using especia::word32;

class Rng_Test : public Unit_Test {
private:

    void test_mt() {
      const word32 seeds[] = {0x123, 0x234, 0x345, 0x456};
      const Mt19937 mt(4, seeds);

      assert_equals(1067595299ul, mt.rand(), "test MT-19937-32 (0)");
      assert_equals(955945823ul, mt.rand(), "test MT-19937-32 (1)");
      assert_equals(477289528ul, mt.rand(), "test MT-19937-32 (2)");
      assert_equals(4107218783ul, mt.rand(), "test MT-19937-32 (3)");
      assert_equals(4228976476ul, mt.rand(), "test MT-19937-32 (4)");
      assert_equals(3344332714ul, mt.rand(), "test MT-19937-32 (5)");
    }

    void test_pcg() {
      const Pcg32 pcg(42ull, 54ull);
      
      assert_equals(0xa15c02b7ul, pcg.rand(), "test PCG-32 (0)");
      assert_equals(0x7b47f409ul, pcg.rand(), "test PCG-32 (1)");
      assert_equals(0xba1d3330ul, pcg.rand(), "test PCG-32 (2)");
      assert_equals(0x83d2f293ul, pcg.rand(), "test PCG-32 (3)");
      assert_equals(0xbfa4784bul, pcg.rand(), "test PCG-32 (4)");
      assert_equals(0xcbed606eul, pcg.rand(), "test PCG-32 (5)");
    }
  
    
    void run_all() override {
        run(this, &Rng_Test::test_mt);
        run(this, &Rng_Test::test_pcg);
    }
};


int main() {
    return Rng_Test().run_testsuite();
}
