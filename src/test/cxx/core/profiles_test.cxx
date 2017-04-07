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

using especia::R_type;
using especia::Equivalent_Width_Calculator;
using especia::Integrator;


class Profiles_Test : public Unit_Test {
private:

    void test_equivalent_width_intergalactic_doppler() {
        using especia::Intergalactic_Doppler;

        const R_type result = calculator.calculate(Intergalactic_Doppler());

        assert_equals(0.698785, result, 1.0E-06, "equivalent width (intergalactic Doppler)");
    }

    void test_equivalent_width_many_multiplet() {
        using especia::Many_Multiplet;

        const R_type result = calculator.calculate(Many_Multiplet());

        assert_equals(0.698785, result, 1.0E-06, "equivalent width (many-multiplet)");
    }

    void test_equivalent_width_intergalactic_voigt() {
        using especia::Intergalactic_Voigt;
        using especia::Pseudo_Voigt;

        const R_type result = calculator.calculate(Intergalactic_Voigt<Pseudo_Voigt>());

        assert_equals(0.928452, result, 1.0E-06, "equivalent width (intergalactic Voigt)");
    }

    void test_equivalent_width_intergalactic_voigt_extended() {
        using especia::Intergalactic_Voigt;
        using especia::Extended_Pseudo_Voigt;

        const R_type result = calculator.calculate(Intergalactic_Voigt<Extended_Pseudo_Voigt>());

        assert_equals(0.928510, result, 1.0E-06, "equivalent width (intergalactic Voigt, extended)");
    }

    void run_all() {
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_doppler);
        run(this, &Profiles_Test::test_equivalent_width_many_multiplet);
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_voigt);
        run(this, &Profiles_Test::test_equivalent_width_intergalactic_voigt_extended);
    }

    Equivalent_Width_Calculator<Integrator<R_type>> calculator;
};


int main() {
    return Profiles_Test().run_testsuite();
}
