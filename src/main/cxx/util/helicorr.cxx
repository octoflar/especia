/// @file helicorr.cxx
/// Utility to apply a heliocentric (or barycentric) velocity correction
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
#include <exception>
#include <stdexcept>

#include "../core/base.h"
#include "../core/dataio.h"
#include "../core/exitcodes.h"

using namespace std;

using especia::N_type;
using especia::R_type;


/**
 * Writes the usage message to an output stream.
 *
 * @param os The ouput stream.
 * @param pname The program name.
 */
void write_usage_message(ostream &os, const string &pname) {
    os << "usage: " << pname << " VELOCITY (m s-1) [SKIP] < ISTREAM > OSTREAM" << endl;
}

/**
 * Utility to apply the heliocentric (or barycentric) velocity correction to
 * spectroscopic data.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The velocity of the observer relative to the heliocenter (or
 * barycenter) of the solar system (m s-1).
 *
 * @c argv[2] The number of lines to skip (optional, default = 0).
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: helicorr VELOCITY (m s-1) [SKIP] < ISTREAM > OSTREAM
 */
int main(int argc, char *argv[]) {
    const string pname(argv[0]);

    if (argc == 1) {
        write_usage_message(cout, pname);
        return 0;
    }
    try {
        if (argc != 2 and argc != 3) {
            throw invalid_argument("Error: an invalid number of arguments was supplied");
        }

        N_type skip = 0;

        if (argc == 3) {
            skip = especia::convert<N_type>(string(argv[2]));
        }

        const R_type v = especia::convert<R_type>(string(argv[1]));

        valarray<R_type> x;
        valarray<R_type> y;
        valarray<R_type> z;

        if (especia::get(cin, x, y, z, skip)) {
            if (v != 0.0) {
                x *= 1.0 + especia::redshift(v);
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
