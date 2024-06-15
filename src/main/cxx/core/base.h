/// @file base.h
/// Base types, mathematical and physical constants and functions.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#ifndef ESPECIA_BASE_H
#define ESPECIA_BASE_H

#include <cinttypes>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

namespace especia
{

/// The type of integer numbers including zero (denoted in maths as set Z).
typedef int integer;

/// The type of natural numbers including zero (denoted in maths as set N).
typedef unsigned natural;

/// The type of real numbers (denoted in maths as set R).
typedef double real;

/// The type of binary numbers with 32 binary digits.
typedef uint32_t word32;

/// The type of binary numbers with 64 binary digits.
typedef uint64_t word64;

/// The class of continuous univariate functions @c f(x) whose derivative
/// exists and is continous.
template <class T = real> class C1
{
public:
  /// The first argument is the abscissa @c x, the second argument returns the
  /// value of @c f(x), and the third argument returns the value of the first
  /// derivative @c f'(x).
  typedef void (&type) (const T &, T &, T &);

private:
  /// The private constructor prevents instantiation.
  C1 () = default;
};

/// Pi. <https://www.wolframalpha.com/input/?i=pi+to+49+digits>
const real pi = real (3.141592653589793238462643383279502884197169399375L);

/// The square root of Pi.
/// <https://www.wolframalpha.com/input/?i=square+root+of+pi+to+49+digits>
const real sqrt_of_pi
    = real (1.772453850905516027298167483341145182797549456123L);

/// The square root of the natural logarithm of 2.
/// <https://www.wolframalpha.com/input/?i=sqrt(ln(2))+to+49+digits>
const real sqrt_of_ln_two
    = real (0.832554611157697756353164644895201047630588852264L);

/// The electric constant (F m-1). *NIST SP 961 (Sept/2015)*
const real electric_constant = real (8.854187817E-12);

/// The electron mass (kg). *NIST SP 961 (Sept/2015)*
const real electron_mass = real (9.10938356E-31);

/// The elementary charge (C). *NIST SP 961 (Sept/2015)*
const real elementary_charge = real (1.6021766208E-19);

/// SI prefix. The spectral resolution of an instrument is expressed in units
/// of this number.
const real kilo = real (1.0E+03);

/// SI prefix.
const real milli = real (1.0E-03);

/// SI prefix. The variation of the fine-structure constant is expressed in
/// units of this number.
const real micro = real (1.0E-06);

/// The speed of light in vacuum (m s-1). *NIST SP 961 (Sept/2015)*
const real speed_of_light = real (299792458.0);

/// Converts a numeric character string into a number.
///
/// @tparam T The number type.
///
/// @param s The string.
/// @return the number.
///
/// @throw invalid_argument when the string cannot be converted into a number
/// of requested type.
template <class T>
T
convert (const std::string &s)
{
  using std::invalid_argument;
  using std::istringstream;

  istringstream iss (s);
  T t;

  if (iss >> t)
    {
      return t;
    }

  throw invalid_argument ("especia::convert(): Error: the expression '" + s
                          + "' cannot be converted into a number");
}

/// Returns the L-2 norm of a vector.
///
/// @tparam T The number type.
///
/// @param[in] n The dimension of the vector.
/// @param[in] x The vector.
/// @return the L-2 norm of the vector.
template <class T = real>
T
norm (natural n, const T x[])
{
  using std::inner_product;
  using std::sqrt;

  return sqrt (inner_product (x, x + n, x, T (0.0)));
}

/// Returns the photon redshift as a function of relative radial velocity
/// between observer and emitter.
///
/// @param[in] v The relative radial velocity between observer and emitter (m
/// s-1).
/// @return the photon redshift.
template <class T = real>
T
redshift (const T &v)
{
  return std::sqrt ((T (1.0) + v / T (speed_of_light))
                    / (T (1.0) - v / T (speed_of_light)))
         - T (1.0);
}

/// Solves the equation f(x) = c by means of Newton's method.
///
/// @tparam T The number type.
///
/// @param[in] f The function.
/// @param[in] c The constant on the right-hand side of the equation.
/// @param[in] x The initial guess of the solution.
/// @param[in] accuracy_goal The accuracy goal.
/// @param[in] max_iteration The maximum number of iterations (optional).
///
/// @return the solution to the equation f(x) = c.
///
/// @throw runtime_error when the accuracy goal was not reached within the
/// prescribed number of iterations.
template <class T = real>
T
solve (typename C1<T>::type &f, T c, T x, T accuracy_goal,
       natural max_iteration = 100)
{
  using std::abs;
  using std::runtime_error;

  T d, y, z;

  for (natural i = 0; i < max_iteration; ++i)
    {
      f (x, y, z);
      d = (y - c) / z;
      x -= d;
      if (abs (d) < accuracy_goal * x)
        {
          return x;
        }
    }

  throw runtime_error (
      "especia::solve() Error: the required accuracy goal was not reached");
}

/// Returns the square of a number.
///
/// @tparam T The number type.
///
/// @param[in] x The number.
/// @return the square of the number.
template <class T = real>
T
sq (const T &x)
{
  return x * x;
}

}

#endif // ESPECIA_BASE_H
