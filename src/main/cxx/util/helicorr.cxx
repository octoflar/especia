/// @file helicorr.cxx
/// Utility to apply a heliocentric (or barycentric) velocity correction
/// Copyright (c) 2016 Ralf Quast
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
#include <cstdlib>

#include "../core/base.h"
#include "../core/dataio.h"

using namespace std;


namespace especia {

    /**
     * Returns the Doppler factor for a given radial velocity.
     *
     * @param[in] v The radial velocity (m s-1).
     * @return the Doppler factor.
     */
    double dopp(const double &v) {
        return std::sqrt((1.0 + v / speed_of_light) / (1.0 - v / speed_of_light));
    }

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
    const char *pname = argv[0];

    int skip = 0;

    if (argc == 3) {
        skip = atoi(argv[2]);
    }
    if (argc == 3 || argc == 2) {
        const double v = atof(argv[1]);

        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        if (especia::get(cin, x, y, z, skip)) {
            if (v != 0.0) {
                x *= especia::dopp(v);
            }
            especia::put(cout, x, y, z);
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " VELOCITY (m s-1) [SKIP] < ISTREAM > OSTREAM" << endl;
        return 1;
    }
}
