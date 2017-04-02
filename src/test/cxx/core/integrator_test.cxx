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


class Integrator_Test {
public:
    void run_testsuite() throw(runtime_error) {
        before_all();
        run_all();
        after_all();
    }

private:
    void before_all() {

    }

    void after_all() {

    }

    void before() {

    }

    void after() {

    }

    template<class T>
    void run(const T &test) {
        before();
        test();
        after();
    }

    void run_all() throw(runtime_error) {
        run(test_integrate_cos);
        run(test_integrate_sin);
        run(test_integrate_sin_sq);
        run(test_integrate_line_profile);
        run(test_integrate_profile__infinite);
    }

    static void test_integrate_cos() {
        using std::cos;
        using especia::Integrator;
        using especia::pi;

        const Assertions assertions;
        const Integrator<double> integrator;

        const double a = 0.0;
        const double b = pi;
        const double result = integrator.integrate([](double x) -> double { return cos(x); }, a, b);

        assertions.assert_equals("test_integrate_cos", 0.0, result, 1.0E-06);
    }

    static void test_integrate_sin() {
        using std::sin;
        using especia::Integrator;
        using especia::pi;

        const Assertions assertions;
        const Integrator<double> integrator;

        const double a = 0.0;
        const double b = pi;
        const double result = integrator.integrate([](double x) -> double { return sin(x); }, a, b);

        assertions.assert_equals("test_integrate_sin", 2.0, result, 1.0E-06);
    }

    static void test_integrate_sin_sq() {
        using std::sin;
        using especia::Integrator;
        using especia::pi;
        using especia::sq;

        const Assertions assertions;
        const Integrator<double> integrator;

        const double a = 0.0;
        const double b = sq(pi);
        const double result = integrator.integrate([](double x) -> double { return sq(sin(x)); }, a, b);

        assertions.assert_equals("test_integrate_sin_sq", 4.740589, result, 1.0E-06);
    }

    static void test_integrate_line_profile() {
        using std::exp;
        using especia::Integrator;
        using especia::pi;
        using especia::sq;

        const Assertions assertions;
        const Integrator<double> integrator;

        const double a = 0.0;
        const double b = 4.0;
        const double result = integrator.integrate([](double x) -> double { return 1.0 - exp(-exp(-sq(x))); }, a, b);

        assertions.assert_equals("test_integrate_profile", 0.642572, result, 1.0E-06);
    }

    static void test_integrate_profile__infinite() {
        using std::exp;
        using std::log;
        using especia::Integrator;
        using especia::pi;
        using especia::sq;

        const Assertions assertions;
        const Integrator<double> integrator;

        const double a = 0.0;
        const double b = 1.0;
        const double result = integrator.integrate([](double u) -> double { return (1.0 - exp(-exp(-sq(log(u))))) / u; }, a, b);

        assertions.assert_equals("test_integrate_profile__infinite", 0.642572, result, 1.0E-06);
    }
};


int main() {
    using std::cout;
    using std::endl;

    Integrator_Test test;

    try {
        test.run_testsuite();
        return 0;
    } catch (runtime_error &e) {
        cout << e.what() << endl;
        return 1;
    }
}
