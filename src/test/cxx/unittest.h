/// @file unittest.h
/// Simple unit testing framework
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#ifndef ESPECIA_UNITTEST_H_H
#define ESPECIA_UNITTEST_H_H

#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

/// The base class to be inherited by all unit-level tests.
class Unit_Test
{
public:
  /// The destructor.
  virtual ~Unit_Test () = default;

  /// Runs the testsuite.
  ///
  /// @return an exit code.
  int
  run_testsuite ()
  {
    using std::endl;
    using std::exception;

    try
      {
        before_all ();
        run_all ();
        after_all ();
        return 0;
      }
    catch (Assertion_Error &e)
      {
        err << e.what () << endl;
        return 1;
      }
    catch (exception &e)
      {
        err << e.what () << endl;
        return 2;
      }
  }

protected:
  /// The type of exception thrown when an assertion fails.
  class Assertion_Error : public std::exception
  {
  public:
    /// The constructor.
    ///
    /// @param what A description of the failed assertion.
    explicit Assertion_Error (std::string what)
        : exception (), what_happened (std::move (what))
    {
    }

    /// The destructor.
    ~Assertion_Error () override = default;

    /// Returns a description of the failed assertion.
    ///
    /// @return the description of the failed assertion.
    const char *
    what () const noexcept override
    {
      return what_happened.c_str ();
    }

  private:
    /// The description of the failed assertion.
    std::string what_happened;
  };

  /// The constructor.
  Unit_Test () = default;

  /// Method called before any test case will be executed.
  virtual void
  before_all ()
  {
  }

  /// Method called after all test cases have been executed.
  virtual void
  after_all ()
  {
  }

  /// Method called before each test case.
  virtual void
  before ()
  {
  }

  /// Method called after each test case.
  virtual void
  after ()
  {
  }

  /// Runs all test cases.
  virtual void run_all () = 0;

  /// Runs a test case.
  ///
  /// @tparam C The test class.
  /// @tparam T The test case type.
  ///
  /// @param c A pointer to the test class.
  /// @param t A pointer to the test case.
  template <class C, class T>
  void
  run (C c, T t)
  {
    before ();
    (c->*t) ();
    after ();
  }

  /// Asserts equality of two values.
  ///
  /// @tparam E The expected type.
  /// @tparam A The actual type.
  ///
  /// @param expected The expected value.
  /// @param actual The actual value.
  /// @param name The assertion name.
  /// @throw an @c Assertion_Error if requested.
  template <class E, class A>
  void
  assert_equals (const E &expected, const A &actual,
                 const std::string &name = "unnamed assertion") const
  {
    using std::cerr;
    using std::endl;
    using std::stringstream;

    if (actual == expected) // NaN safe
      {
        handle_passed (name);
      }
    else
      {
        stringstream what;
        what << "Failed: " << name << "\n";
        what << "    Expected result: " << expected << "\n";
        what << "    Actual   result: " << actual;
        handle_failed (what.str ());
      }
  }

  /// Asserts equality of two values with some (absolute) tolerance.
  ///
  /// @tparam E The expected type.
  /// @tparam A The actual type.
  ///
  /// @param expected The expected value.
  /// @param actual The actual value.
  /// @param tolerance The (absolute) tolerance.
  /// @param name The assertion name.
  /// @throw an @c Assertion_Error if requested.
  template <class E, class A>
  void
  assert_equals (const E &expected, const A &actual, const E &tolerance,
                 const std::string &name = "unnamed assertion") const
  {
    using std::abs;
    using std::stringstream;

    if (abs (actual - expected) <= tolerance) // NaN safe
      {
        handle_passed (name);
      }
    else
      {
        stringstream what;
        what << "Failed: " << name << "\n";
        what << "    Expected result: " << expected << "\n";
        what << "    Actual   result: " << actual << "\n";
        what << "    Expected tolerance: " << tolerance;
        handle_failed (what.str ());
      }
  }

  /// Asserts @c false.
  ///
  /// @param actual The actual value.
  /// @param name The assertion name.
  /// @throw an @c Assertion_Error if requested.
  void
  assert_false (const bool actual,
                const std::string &name = "unnamed assertion") const
  {
    assert_equals (false, actual, name);
  }

  /// Asserts @c true.
  ///
  /// @param actual The actual value.
  /// @param name The assertion name.
  /// @throw an @c Assertion_Error if requested.
  void
  assert_true (const bool actual,
               const std::string &name = "unnamed assertion") const
  {
    assert_equals (true, actual, name);
  }

private:
  /// Handles a failed assertion.
  ///
  /// @param what A description of the failed assertion.
  /// @throw an @c Assertion_Error if requested.
  void
  handle_failed (const std::string &what) const
  {
    using std::endl;

    if (message_on_failed)
      {
        err << what << endl;
      }
    if (throw_on_failed)
      {
        throw Assertion_Error (what); // NOLINT
      }
  }

  /// Handles a passed assertion.
  ///
  /// @param name The name of the passed assertion.
  void
  handle_passed (const std::string &name) const
  {
    using std::endl;

    if (message_on_passed)
      {
        out << "Passed: " << name << endl;
      }
  }

  /// Issue a message when an assertion failed?
  bool message_on_failed = false;

  /// Issue a message when an assertion passed?
  bool message_on_passed = true;

  /// Throw an exception when an assertion failed?
  bool throw_on_failed = true;

  /// The output stream for error messages.
  std::ostream &err = std::cerr;

  /// The output stream for other messages.
  std::ostream &out = std::cout;
};

#endif // ESPECIA_UNITTEST_H_H
