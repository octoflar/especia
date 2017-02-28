// Especia
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
#include <exception>
#include <iomanip>
#include <iostream>

#include "config.h"
#include "core/decompose.h"
#include "core/model.h"
#include "core/mtwister.h"
#include "core/randev.h"

const char usemsg[] = "usage: ";
const char parmsg[] = "SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE < ISTREAM > OSTREAM";

/**
 * Flavor of Especia to analyse intergalactic Lyman-alpha lines.
 *
 * @param argc The number of command line arguments supplied.
 * @param argv The command line arguments:
 * @parblock
 * @c argv[0] The program name
 *
 * @c argv[1] The random seed
 *
 * @c argv[2] The parent number
 *
 * @c argv[3] The population size
 *
 * @c argv[4] The initial global step size
 *
 * @c argv[5] The accuracy goal
 *
 * @c argv[6] The stop generation number
 *
 * @c argv[7] The trace interval
 * @endparblock
 * @return an exit code.
 *
 * @remark Usage: especiv SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE < ISTREAM > OSTREAM
 */
int main(int argc, char *argv[]) {
    using namespace especia;
    using std::cin;
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::exception;

    const char *pname = argv[0];

    int exit_code = 0;

    if (argc == 8) {
        const unsigned long seed = (unsigned long) atol(argv[1]);
        const unsigned parent_number = (unsigned) atoi(argv[2]);
        const unsigned population_size = (unsigned) atoi(argv[3]);
        const double step_size = atof(argv[4]);
        const double accuracy_goal = atof(argv[5]);
        const unsigned long stop_generation = (unsigned long) atol(argv[6]);
        const unsigned trace = (unsigned) atoi(argv[7]);

        cout << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
        cout << "<html>\n";
        cout << "<!--\n";
        cout << "<command>\n";
        for (int i = 0; i < argc; ++i) {
            cout << " " << argv[i];
        }
        cout << endl;
        cout << "</command>\n";
        cout << "-->\n";
        cout << "</html>\n";

        Model<G_Voigt<Pseudo_Voigt>> model;
        model.get(cin, cout);

        if (cin.eof() and !cin.fail()) {
            try {
                Normal_Deviate<MT19937> normal_deviate(seed);
                Decompose decompose;

                if (model.optimize(parent_number,
                                   population_size,
                                   step_size,
                                   accuracy_goal,
                                   stop_generation,
                                   trace,
                                   normal_deviate, decompose, cout)) {
                    exit_code = 0;
                } else {
                    exit_code = 2;
                }
                model.put(cout);
            } catch (exception &e) {
                cerr << e.what() << endl;

                exit_code = 3;
            }
        }
    } else {
        cout << PROJECT_LONG_NAME << " " << DOI << endl;
        cout << usemsg << pname << ": " << parmsg << endl;

        exit_code = 1;
    }

    return exit_code;
}
