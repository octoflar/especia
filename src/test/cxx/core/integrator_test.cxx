/// @file integrator_test.cxx
/// Unit tests
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include "../../../main/cxx/core/base.h"
#include "../../../main/cxx/core/integrator.h"
#include "../unittest.h"

using especia::Integrator;
using especia::real;

class Integrator_Test : public Unit_Test
{
public:
  /// Tests the default integrator.
  Integrator_Test () : integrator () {}

  /// Tests a custom integrator.
  Integrator_Test (Integrator<>::Formula p, Integrator<>::Formula q)
      : integrator (p, q)
  {
  }

  ~Integrator_Test () = default;

private:
  void
  test_integrate_constant ()
  {
    const real result
        = integrator.integrate ([] (real x) -> real { return real (1.0); },
                                real (0.0), real (1.0), real (1.0E-12));

    assert_equals (real (1.0), result, real (0.5E-12), "integrate constant");
  }

  void
  test_integrate_identity ()
  {
    const real result
        = integrator.integrate ([] (real x) -> real { return x; }, real (0.0),
                                real (1.0), real (1.0E-12));

    assert_equals (real (0.5), result, real (0.5E-12), "integrate identity");
  }

  void
  test_integrate_parabola ()
  {
    const real result
        = integrator.integrate ([] (real x) -> real { return x * x; },
                                real (0.0), real (1.0), real (1.0E-12));

    assert_equals (real (1.0 / 3.0), result, real (0.5E-12),
                   "integrate parabola");
  }

  void
  test_integrate_cos ()
  {
    using especia::pi;
    using std::cos;

    const real result = integrator.integrate (
        [] (real x) -> real { return cos (x); }, real (0.0), pi);

    assert_equals (real (0.0), result, real (0.5E-06), "integrate cosine");
  }

  void
  test_integrate_sin ()
  {
    using especia::pi;
    using std::sin;

    const real result = integrator.integrate (
        [] (real x) -> real { return sin (x); }, real (0.0), pi);

    assert_equals (real (2.0), result, real (0.5E-06), "integrate sine");
  }

  void
  test_integrate_sin_sq ()
  {
    using especia::pi;
    using especia::sq;
    using std::sin;

    const real result
        = integrator.integrate ([] (real x) -> real { return sq (sin (x)); },
                                real (0.0), real (2.0) * pi);

    assert_equals (pi, result, real (0.5E-06), "integrate sine squared");
  }

  void
  test_integrate_absorption ()
  {
    using especia::sq;
    using std::exp;

    const real result = integrator.integrate (
        [] (real x) -> real { return real (1.0) - exp (-exp (-sq (x))); },
        real (0.0), real (4.0));
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+4%7D%5D>
    assert_equals (real (0.642572), result, real (0.5E-06),
                   "integrate absorption");
  }

  void
  test_integrate_absorption_positive_infinite ()
  {
    using especia::sq;
    using std::exp;

    const real result = integrator.integrate_positive_infinite (
        [] (real x) -> real { return real (1.0) - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+0,+Infinity%7D%5D>
    assert_equals (real (0.642572), result, real (0.5E-06),
                   "integrate absorption (positive-infinite)");
  }

  void
  test_integrate_absorption_negative_infinite ()
  {
    using especia::sq;
    using std::exp;

    const real result = integrator.integrate_negative_infinite (
        [] (real x) -> real { return real (1.0) - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+0%7D%5D>
    assert_equals (real (0.642572), result, real (0.5E-06),
                   "integrate absorption (negative-infinite");
  }

  void
  test_integrate_absorption_infinite ()
  {
    using especia::sq;
    using std::exp;

    const real result = integrator.integrate_infinite (
        [] (real x) -> real { return real (1.0) - exp (-exp (-sq (x))); });
    // <https://www.wolframalpha.com/input/?i=integrate%5B1-Exp%5B-Exp%5B-x%5E2%5D%5D,%7Bx,+-Infinity,+Infinity%7D%5D>
    assert_equals (real (1.285145), result, real (0.5E-06),
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

  const Integrator<real> integrator;
};

int
main ()
{
  const auto Q13 = Integrator<>::Formula::Q13;
  const auto Q19 = Integrator<>::Formula::Q19;
  const auto Q27 = Integrator<>::Formula::Q27;
  const auto Q41 = Integrator<>::Formula::Q41;

  return Integrator_Test ().run_testsuite () // test the default first
         or Integrator_Test (Q13, Q19).run_testsuite ()
         or Integrator_Test (Q13, Q27).run_testsuite ()
         or Integrator_Test (Q13, Q41).run_testsuite ()
         or Integrator_Test (Q19, Q27).run_testsuite ()
         or Integrator_Test (Q19, Q41).run_testsuite ()
         or Integrator_Test (Q27, Q41).run_testsuite ();
}
