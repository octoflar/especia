/// @file profiles_test.cxx
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
#include "../../../main/cxx/core/base.h"
#include "../../../main/cxx/core/integrator.h"
#include "../../../main/cxx/core/profiles.h"
#include "../unittest.h"

using especia::real;
using especia::Equivalent_Width_Calculator;
using especia::Integrator;


class Profiles_Test : public Unit_Test {
private:

    void test_equivalent_width_intergalactic_doppler() {
        using especia::Intergalactic_Doppler;

        const real result = calculator.calculate(Intergalactic_Doppler()) / real(1000.0);

        // <https://www.wolframalpha.com/input/?i=NIntegrate%5B1-Exp%5B-1%2F(0.5+Sqrt%5BPi%5D)+Exp%5B-(x%2F0.5)%5E2%5D%5D,+%7Bx,+-Infinity,+Infinity%7D%5D>
        assert_equals(real(0.698785), result, real(1.0E-06), "equivalent width (intergalactic Doppler)");
    }

    void test_equivalent_width_many_multiplet() {
        using especia::Many_Multiplet;

        const real result = calculator.calculate(Many_Multiplet()) / real(1000.0);

        // <https://www.wolframalpha.com/input/?i=NIntegrate%5B1-Exp%5B-1%2F(0.5+Sqrt%5BPi%5D)+Exp%5B-(x%2F0.5)%5E2%5D%5D,+%7Bx,+-Infinity,+Infinity%7D%5D>
        assert_equals(real(0.698785), result, real(1.0E-06), "equivalent width (many-multiplet)");
    }

    void test_equivalent_width_intergalactic_voigt() {
        using especia::Intergalactic_Voigt;
        using especia::Pseudo_Voigt;

        const real result = calculator.calculate(Intergalactic_Voigt<Pseudo_Voigt>()) / real(1000.0);

        // <https://www.wolframalpha.com/input/?i=NIntegrate%5B1+-+Exp%5B-PDF%5BVoigtDistribution%5B0.5,+0.5%2FSqrt%5B2%5D%5D,+x%5D%5D,+%7Bx,+-Infinity,+Infinity%7D%5D>
        assert_equals(real(0.881143), result, real(3.5E-03), "equivalent width (intergalactic Voigt)");
    }

    void test_equivalent_width_intergalactic_voigt_extended() {
        using especia::Intergalactic_Voigt;
        using especia::Extended_Pseudo_Voigt;

        const real result = calculator.calculate(Intergalactic_Voigt<Extended_Pseudo_Voigt>()) / real(1000.0);

        // <https://www.wolframalpha.com/input/?i=NIntegrate%5B1+-+Exp%5B-PDF%5BVoigtDistribution%5B0.5,+0.5%2FSqrt%5B2%5D%5D,+x%5D%5D,+%7Bx,+-Infinity,+Infinity%7D%5D>
        assert_equals(real(0.881143), result, real(4.5E-03), "equivalent width (intergalactic Voigt, extended)");
    }

    void test_maximum_pseudo_voigt() {
        using especia::Pseudo_Voigt;

        // <https://www.wolframalpha.com/input/?i=PDF%5BVoigtDistribution%5B0.5,+0.5%2FSqrt%5B2%5D%5D,+0%5D>
        assert_equals(real(0.482476), Pseudo_Voigt(0.5, 0.5)(0.0), real(1.0E-03),
                      "Voigt function maximum (pseudo-Voigt approximation)");

        // <https://www.wolframalpha.com/input/?i=PDF%5BVoigtDistribution%5B1.0,+1.0%2FSqrt%5B2%5D%5D,+0%5D>
        assert_equals(real(0.241238), Pseudo_Voigt(1.0, 1.0)(0.0), real(1.0E-03),
                      "Voigt function maximum (pseudo-Voigt approximation)");
    }

    void test_maximum_pseudo_voigt_extended() {
        using especia::Extended_Pseudo_Voigt;

        // <https://www.wolframalpha.com/input/?i=PDF%5BVoigtDistribution%5B0.5,+0.5%2FSqrt%5B2%5D%5D,+0%5D>
        assert_equals(real(0.482476), Extended_Pseudo_Voigt(0.5, 0.5)(0.0), real(0.5E-03),
                      "Voigt function maximum (extended pseudo-Voigt approximation)");

        // <https://www.wolframalpha.com/input/?i=PDF%5BVoigtDistribution%5B1.0,+1.0%2FSqrt%5B2%5D%5D,+0%5D>
        assert_equals(real(0.241238), Extended_Pseudo_Voigt(1.0, 1.0)(0.0), real(0.5E-03),
                      "Voigt function maximum (extended pseudo-Voigt approximation)");
    }

    void run_all() override {
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_doppler);
        run(this, &Profiles_Test::test_equivalent_width_many_multiplet);
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_voigt);
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_voigt_extended);
        run(this, &Profiles_Test::test_maximum_pseudo_voigt);
        run(this, &Profiles_Test::test_maximum_pseudo_voigt_extended);
    }

    Equivalent_Width_Calculator<Integrator<real>> calculator;
};


int main() {
    return Profiles_Test().run_testsuite();
}
