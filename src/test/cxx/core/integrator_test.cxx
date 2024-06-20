/// @file integrator_test.cxx
/// Unit tests
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include "../../../main/cxx/core/base.h"
#include "../../../main/cxx/core/integrator.h"
#include "../unittest.h"

using especia::Integrator;

class Integrator_Test : public Unit_Test
{
private:
  void
  test_integrate_constant ()
  {
    const double result = integrator.integrate (
        [] (double x) -> double { return 1.0; }, 0.0, 1.0, 1.0E-12);

    assert_equals (1.0, result, 1.0E-12, "integrate constant");
  }

  void
  test_integrate_identity ()
  {
    const double result = integrator.integrate (
        [] (double x) -> double { return x; }, 0.0, 1.0, 1.0E-12);

    assert_equals (0.5, result, 1.0E-12, "integrate identity");
  }

  void
  test_integrate_parabola ()
  {
    const double result = integrator.integrate (
        [] (double x) -> double { return x * x; }, 0.0, 1.0, 1.0E-12);

    assert_equals (1.0 / 3.0, result, 1.0E-12, "integrate parabola");
  }

  void
  test_integrate_cos ()
  {
    using especia::pi;
    using std::cos;

    const double result = integrator.integrate (
        [] (double x) -> double { return cos (x); }, 0.0, pi);

    assert_equals (0.0, result, 0.5E-06, "integrate cosine");
  }

  void
  test_integrate_sin ()
  {
    using especia::pi;
    using std::sin;

    const double result = integrator.integrate (
        [] (double x) -> double { return sin (x); }, 0.0, pi);

    assert_equals (2.0, result, 0.5E-06, "integrate sine");
  }

  void
  test_integrate_sin_sq ()
  {
    using especia::pi;
    using especia::sq;
    using std::sin;

    const double result = integrator.integrate (
        [] (double x) -> double { return sq (sin (x)); }, 0.0, 2.0 * pi);

    assert_equals (pi, result, 0.5E-06, "integrate sine squared");
  }

  void
  test_integrate_absorption ()
  {
    using especia::sq;
    using std::exp;

    const double result = integrator.integrate (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); }, 0.0,
        4.0);
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+4%7D%5D>
    assert_equals (0.642572, result, 0.5E-06, "integrate absorption");
  }

  void
  test_integrate_absorption_positive_infinite ()
  {
    using especia::sq;
    using std::exp;

    double result;

    result = q13.integrate_positive_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+Infinity%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (positive-infinite, Q13-Q19)");

    result = q19.integrate_positive_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+Infinity%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (positive-infinite, Q19-Q27)");

    result = q27.integrate_positive_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+Infinity%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (positive-infinite, Q27-Q41)");
  }

  void
  test_integrate_absorption_negative_infinite ()
  {
    using especia::sq;
    using std::exp;

    double result;

    result = q13.integrate_negative_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+0%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (negative-infinite, Q13-Q19)");

    result = q19.integrate_negative_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+0%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (negative-infinite, Q19-Q27)");

    result = q27.integrate_negative_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+0%7D%5D>
    assert_equals (0.642572, result, 0.5E-06,
                   "integrate absorption (negative-infinite, Q27-Q41)");
  }

  void
  test_integrate_absorption_infinite ()
  {
    using especia::sq;
    using std::exp;

    const double result = integrator.integrate_infinite (
        [] (double x) -> double { return 1.0 - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+Infinity%7D%5D>
    assert_equals (1.285145, result, 0.5E-06,
                   "integrate absorption (infinite)");
  }

  void
  run_all () override
  {
    run (this, &Integrator_Test::test_integrate_constant);
    run (this, &Integrator_Test::test_integrate_identity);
    run (this, &Integrator_Test::test_integrate_parabola);
    run (this, &Integrator_Test::test_integrate_cos);
    run (this, &Integrator_Test::test_integrate_sin);
    run (this, &Integrator_Test::test_integrate_sin_sq);
    run (this, &Integrator_Test::test_integrate_absorption);
    run (this, &Integrator_Test::test_integrate_absorption_positive_infinite);
    run (this, &Integrator_Test::test_integrate_absorption_negative_infinite);
    run (this, &Integrator_Test::test_integrate_absorption_infinite);
  }

  /// The default integrator
  const Integrator<double> integrator;

  const Integrator<double> q13
      = Integrator<> (Integrator<>::Formula::Q13, Integrator<>::Formula::Q19);
  const Integrator<double> q19
      = Integrator<> (Integrator<>::Formula::Q19, Integrator<>::Formula::Q27);
  const Integrator<double> q27
      = Integrator<> (Integrator<>::Formula::Q27, Integrator<>::Formula::Q41);
};

int
main ()
{
  return Integrator_Test ().run_testsuite ();
}
