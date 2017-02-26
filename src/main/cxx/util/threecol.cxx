// Utility: merge separate spectral flux and uncertainty data to three-column format
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
#include <fstream>

#include "dataio.h"

using namespace std;


/**
 * Utility to merge separate spectral flux and uncertainty data to
 * three-column format.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The path name of the flux data file.
 *
 * @c argv[2] The path name of the flux uncertainty data file.
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: threecol FLUX UNCERTAINTY > OSTREAM
 */
int main(int argc, char *argv[]) {
    const char *pname = argv[0];

    if (argc == 3) {
        valarray<double> x;
        valarray<double> y;
        valarray<double> z;

        ifstream fxy(argv[1]);
        ifstream fxz(argv[2]);

        especia::get(fxy, x, y);
        especia::get(fxz, x, z);

        fxy.close();
        fxz.close();

        if (fxy and fxz and y.size() == z.size()) {
            especia::put(cout, x, y, z);
        } else {
            cerr << pname << ": input failure" << endl;
            return 2;
        }

        return 0;
    } else {
        cout << "usage: " << pname << " FLUX UNCERTAINTY > OSTREAM" << endl;
        return 1;
    }
}
