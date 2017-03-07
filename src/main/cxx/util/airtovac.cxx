/// @file airtovac.cxx
/// Utility to convert photon wavelength from air to vacuum
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
#include <cmath>
#include <exception>
#include <stdexcept>

#include "../core/base.h"
#include "../core/dataio.h"
#include "../core/exitcodes.h"

using namespace std;

using especia::N_type;
using especia::R_type;


/**
 * Solves the equation f(x) = c by means of Newton's method.
 *
 * @param[in] f A differentiable function, which takes as parameters:
 * @parblock
 * @param[in] x The abscissa value.
 * @param[out] y The result y = f(x).
 * @param[out] z The derivative of @c y with respect to @c x.
 * @endparblock
 * @param[in] c The constant on the right-hand side of the equation.
 * @param[in] x The initial guess of the solution.
 * @param[in] accuracy_goal The accuracy goal (optional).
 *
 * @return the solution to the equation f(x) = c.
 */
R_type solve(void f(const R_type &x, R_type &y, R_type &z), R_type c, R_type x,
             R_type accuracy_goal = 1.0E-8) throw(std::runtime_error) {
    using std::abs;
    using std::runtime_error;

    R_type d, y, z;

    for (N_type i = 0; i < 100; ++i) {
        f(x, y, z);
        d = (y - c) / z;
        x -= d;
        if (abs(d) < accuracy_goal * x) {
            return x;
        }
    }

    throw runtime_error("Error: the required accuracy goal was not reached");
}

/**
 * Utility to convert photon wavelength (Angstrom) in spectroscopic data from air
 * to vacuum.
 *
 * Further reading:
 *
 * B. Edlen (1966).
 *   *The refractive index of air.*
 *   Metrologia, 2, 2, 71-80.
 *   http://dx.doi.org/10.1088/0026-1394/2/2/002
 *
 * B. Edlen (1953).
 *   *The dispersion of standard air.*
 *   Journal of the Optical Society of America, 43, 5, 339.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The number of lines to skip (optional, default = 0).
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: airtovac [SKIP] < ISTREAM > OSTREAM
 */
int main(int argc, char *argv[]) {
    const string pname(argv[0]);

    try {
        if (argc != 1 and argc != 2) {
            throw invalid_argument("Error: an invalid number of arguments was supplied");
        }

        N_type skip = 0;

        if (argc == 2) {
            skip = especia::convert<N_type>(string(argv[1]));
        }

        valarray<R_type> x;
        valarray<R_type> y;
        valarray<R_type> z;

        if (especia::get(cin, x, y, z, skip)) {
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] = 10.0 / solve(especia::edlen66, 10.0 / x[i], 10.0 / x[i]);
            }
            especia::put(cout, x, y, z);
        } else {
            throw runtime_error("Error: an input error occurred");
        }
        return 0;
    } catch (logic_error &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::LOGICAL_ERROR;
    } catch (runtime_error &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::RUNTIME_ERROR;
    } catch (exception &e) {
        cerr << e.what() << endl;
        return especia::Exit_Codes::UNSPECIFIC_EXCEPTION;
    }
}
