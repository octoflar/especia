// Utility: convert wavelength in spectroscopic data from air to vacuum
// Copyright (c) 2016 Ralf Quast
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#include <cstdlib>
#include <exception>

#include "../base.h"
#include "../dataio.h"

using namespace std;

using especia::edlen;
using especia::get;
using especia::put;


double
solve(void f(const double &, double &, double &), double c, double x, double accuracy_goal) throw(runtime_error) {
    double d, y, z;

    for (unsigned int i = 0; i < 100; ++i) {
        f(x, y, z);
        d = (y - c) / z;
        x -= d;
        if (abs(d) < accuracy_goal * x) {
            return x;
        }
    }

    throw runtime_error("solve(): Error: failed to reach accuracy goal");
}

/**
 * Utility to convert wavelength in spectroscopic data from air
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
 *
 * @c argv[2] The accuracy goal (optional, dafault = 1.0E-8).
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: airtovac [SKIP] [ACCURACY] < ISTREAM > OSTREAM
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    int skip = 0;
    double accuracy_goal = 1.0e-08;

    if (argc == 3) {
        accuracy_goal = atof(argv[2]);
    }
    if (argc == 3 || argc == 2) {
        skip = atoi(argv[1]);
    }
    if (argc == 3 || argc == 2 || argc == 1) {
        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        if (get(cin, x, y, z, skip)) {
            try {
                for (size_t i = 0; i < x.size(); ++i) {
                    x[i] = 10.0 / solve(edlen, 10.0 / x[i], 10.0 / x[i], accuracy_goal);
                }
            }
            catch (exception &e) {
                cerr << pname << ": " << e.what() << endl;
                return 3;
            }
            put(cout, x, y, z);
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " [SKIP] [ACCURACY] < ISTREAM > OSTREAM" << endl;
        return 1;
    }
}
