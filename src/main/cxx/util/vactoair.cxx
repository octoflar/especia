// Utility: convert photon wavelength in spectroscopic data from vacuum to air
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

#include "dataio.h"

using namespace std;


namespace especia {

    /**
     * Returns the square of a number.
     *
     * @tparam Number The number type.
     *
     * @param x[in] The number.
     * @return the square of the number.
     */
    template<class Number>
    inline
    Number sqr(const Number &x) {
        return (x == Number(0)) ? Number(0) : x * x;
    }

    /**
     * Used to convert photon wavelength in vacuum to photon wavelength in air.
     *
     * Further reading:
     *
     * B. Edlén (1966).
     *   *The refractive index of air.*
     *   Metrologia, 2, 2, 71-80.
     *   http://dx.doi.org/10.1088/0026-1394/2/2/002
     *
     * @param x[in] The wavenumber in vacuum (nm-1).
     * @return the wavenumber in air (nm-1).
     *
     * @attention The function uses wavenumber (nm-1) := 10.0 / wavelength (Angstrom) as input and output.
     */
    inline
    double edlen_1966(const double &x) {
        return (1.0000834213 + 1.5997E-10 / (0.0000389 - x * x) + 2.406030E-08 / (0.000130 - x * x)) * x;
    }

}

/**
 * Utility to convert photon wavelength (Angstrom) in spectroscopic data from
 * vacuum to air.
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
 * @remark Usage: vactoair [SKIP] < ISTREAM > OSTREAM
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    int skip = 0;

    if (argc == 2) {
        skip = atoi(argv[1]);
    }
    if (argc == 2 || argc == 1) {
        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        if (especia::get(cin, x, y, z, skip)) {
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] = 10.0 / especia::edlen_1966(10.0 / x[i]);
            }
            especia::put(cout, x, y, z);
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " [SKIP] < ISTREAM > OSTREAM" << endl;
        return 1;
    }
}
