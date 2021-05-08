/// @file integrator_test.cxx
/// Unit tests
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include "../../../main/cxx/core/base.h"
#include "../../../main/cxx/core/integrator.h"
#include "../unittest.h"

using especia::Integrator;


class Integrator_Test : public Unit_Test {
private:

    void test_integrate_cos() {
        using std::cos;
        using especia::pi;

        const double result = integrator.integrate([](double x) -> double { return cos(x); }, 0.0, pi);

        assert_equals(0.0, result, 0.5E-06, "integrate cosine");
    }

    void test_integrate_sin() {
        using std::sin;
        using especia::pi;

        const double result = integrator.integrate([](double x) -> double { return sin(x); }, 0.0, pi);

        assert_equals(2.0, result, 0.5E-06, "integrate sine");
    }

    void test_integrate_sin_sq() {
        using std::sin;
        using especia::pi;
        using especia::sq;

        const double result = integrator.integrate([](double x) -> double { return sq(sin(x)); }, 0.0, 2.0 * pi);

        assert_equals(pi, result, 0.5E-06, "integrate sine squared");
    }

    void test_integrate_absorption() {
        using std::exp;
        using especia::sq;

        const double result = integrator.integrate(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); }, 0.0, 4.0);

        // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+4%7D%5D>
        assert_equals(0.642572, result, 0.5E-06, "integrate absorption");
    }

    void test_integrate_absorption_positive_infinite() {
        using std::exp;
        using especia::sq;

        const double result = integrator.integrate_positive_infinite(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); });

        // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+Infinity%7D%5D>
        assert_equals(0.642572, result, 0.5E-06, "integrate absorption (positive-infinite)");
    }

    void test_integrate_absorption_negative_infinite() {
        using std::exp;
        using especia::sq;

        const double result = integrator.integrate_negative_infinite(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); });

        // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+0%7D%5D>
        assert_equals(0.642572, result, 0.5E-06, "integrate absorption (negative-infinite)");
    }

    void test_integrate_absorption_infinite() {
        using std::exp;
        using especia::sq;

        const double result = integrator.integrate_infinite(
                [](double x) -> double { return 1.0 - exp(-exp(-sq(x))); });

        // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+Infinity%7D%5D>
        assert_equals(1.285145, result, 0.5E-06, "integrate absorption (infinite)");
    }

    void run_all() override {
        run(this, &Integrator_Test::test_integrate_cos);
        run(this, &Integrator_Test::test_integrate_sin);
        run(this, &Integrator_Test::test_integrate_sin_sq);
        run(this, &Integrator_Test::test_integrate_absorption);
        run(this, &Integrator_Test::test_integrate_absorption_positive_infinite);
        run(this, &Integrator_Test::test_integrate_absorption_negative_infinite);
        run(this, &Integrator_Test::test_integrate_absorption_infinite);
    }

    Integrator<double> integrator;
};


int main() {
    return Integrator_Test().run_testsuite();
}
