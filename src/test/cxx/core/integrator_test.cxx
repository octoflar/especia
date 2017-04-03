/// @file integrator_test.cxx
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
#include "../unittest.h"

using especia::Integrator;


class Integrator_Test : public Unit_Test {
private:

    void test_integrate_cos() const {
        using std::cos;
        using especia::pi;

        const double result = integrator.integrate([](double x) -> double { return cos(x); }, 0.0, pi);

        assert_equals(0.0, result, 1.0E-06, "integrate cosine");
    }

    void test_integrate_sin() {
        using std::sin;
        using especia::pi;

        const double result = integrator.integrate([](double x) -> double { return sin(x); }, 0.0, pi);

        assert_equals(2.0, result, 1.0E-06, "integrate sine");
    }

    void test_integrate_sin_sq() {
        using std::sin;
        using especia::pi;
        using especia::sq;

        const double result = integrator.integrate([](double x) -> double { return sq(sin(x)); }, 0.0, sq(pi));

        assert_equals(4.740589, result, 1.0E-06, "integrate sine squared");
    }

    void test_integrate_optical_depth() {
        using std::exp;
        using especia::sq;

        const double result = integrator.integrate(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); }, 0.0, 4.0);

        assert_equals(0.642572, result, 1.0E-06, "integrate optical depth");
    }

    void test_integrate_optical_depth_semi_infinite() {
        using std::exp;
        using std::log;
        using especia::sq;

        const double result = integrator.integrate_semi_infinite(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); });

        assert_equals(0.642572, result, 1.0E-06, "integrate optical depth (semi-infinite)");
    }

    void run_all() {
        run(this, &Integrator_Test::test_integrate_cos);
        run(this, &Integrator_Test::test_integrate_sin);
        run(this, &Integrator_Test::test_integrate_sin_sq);
        run(this, &Integrator_Test::test_integrate_optical_depth);
        run(this, &Integrator_Test::test_integrate_optical_depth_semi_infinite);
    }

    Integrator<double> integrator;
};


int main() {
    using std::cerr;
    using std::endl;

    Integrator_Test test;

    try {
        test.run_testsuite();
        return 0;
    } catch (std::exception &e) {
        cerr << e.what() << endl;
        return 1;
    }
}
