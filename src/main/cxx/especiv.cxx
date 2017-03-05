/// @file especiv.cxx
/// Especia for intergalactic metal and damped H I, He I, II lines.
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
#include <exception>
#include <iostream>

#include "core/model.h"
#include "core/profiles.h"
#include "runner/runner.h"

using namespace std;


/**
 * Flavor of Especia to analyse intergalactic metal and damped H I, He I, II lines.
 *
 * @param argc The number of command line arguments.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name.
 *
 * @c argv[1] The random seed.
 *
 * @c argv[2] The parent number.
 *
 * @c argv[3] The population size.
 *
 * @c argv[4] The initial global step size.
 *
 * @c argv[5] The accuracy goal.
 *
 * @c argv[6] The stop generation number.
 *
 * @c argv[7] The trace interval.
 * @endparblock
 *
 * @return an exit code: 0 = OK, 1 = model is not optimized, 2-4 = input error, 5 = unspecified error
 *
 * @remark Usage: especiv SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE < ISTREAM > OSTREAM
 *
 * @remark A usage message is witten to standard output, if no command line arguments (excluding the
 * program name) are supplied. In this case the returned exit code is zero.
 */
int main(int argc, char *argv[]) {
    especia::Model<especia::Intergalactic_Voigt<especia::Pseudo_Voigt>> model;

    try {
        return especia::Runner(argc, argv).run(model);
    } catch (exception &e) {
        cerr << e.what() << endl;

        return 5;
    }
}
