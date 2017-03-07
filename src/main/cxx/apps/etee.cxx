/// @file etee.cxx
/// Utility to merge separated spectral flux and uncertainty data
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
#include <fstream>

#include "../core/base.h"
#include "../core/dataio.h"

using namespace std;

using especia::Nint_t;
using especia::Real_t;


void write_usage_message(ostream &os, const string &pname) {
    os << "usage: " << pname << " FLUX UNCERTAINTY [SKIP] > OSTREAM" << endl;
}

/**
 * Utility to merge separated spectral flux and uncertainty data columns.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The path name of the flux data file.
 *
 * @c argv[2] The path name of the flux uncertainty data file.
 *
 * @c argv[3] The number of lines to skip (optional, default = 0).
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: etee FLUX UNCERTAINTY [SKIP] > OSTREAM
 */
int main(int argc, char *argv[]) {
    const string pname(argv[0]);

    if (argc == 1) {
        write_usage_message(cout, pname);
        return 0;
    }
    try {
        if (argc != 3 and argc != 4) {
            throw invalid_argument("Error: an invalid number of arguments was supplied");
        }

        Nint_t skip = 0;

        if (argc == 4) {
            skip = especia::convert<Nint_t>(string(argv[2]));
        }

        valarray<Real_t> x;
        valarray<Real_t> y;
        valarray<Real_t> z;

        ifstream fxy(argv[1]);
        ifstream fxz(argv[2]);

        especia::get(fxy, x, y, skip);
        especia::get(fxz, x, z, skip);

        fxy.close();
        fxz.close();

        if (fxy and fxz and y.size() == z.size()) {
            especia::put(cout, x, y, z);
        } else {
            throw runtime_error("Error: an input error occurred");
        }
        return 0;
    } catch (invalid_argument &e) {
        cerr << e.what() << endl;
        write_usage_message(cout, pname);
        return 10;
    } catch (runtime_error &e) {
        cerr << e.what() << endl;
        return 20;
    } catch (exception &e) {
        cerr << e.what() << endl;
        return 30;
    }
}
