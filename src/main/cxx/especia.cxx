// Especia: many-multiplets version to infer the variation of the fine-structure
// constant
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

#define RQ_MANY_MULTIPLET_ANALYSIS 1

#include "model.h"

#undef RQ_MANY_MULTIPLET_ANALYSIS

#include "mtwister.h"
#include "randev.h"
#include "symeig.h"

const char usemsg[] = "usage: ";
const char parmsg[] = "SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE < ISTREAM > OSTREAM";

int main(int argc, char *argv[]) {
    using namespace RQ;
    using std::cin;
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::exception;

    const char *pname = argv[0];

    if (argc == 8) {
        const long seed = atol(argv[1]);
        const int parent_number = atoi(argv[2]);
        const int population_size = atoi(argv[3]);
        const double step_size = atof(argv[4]);
        const double accuracy_goal = atof(argv[5]);
        const long stop_generation = atol(argv[6]);
        const int trace = atoi(argv[7]);

        cout << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
        cout << "<html>\n";
        cout << "<!--\n";
        cout << "<command>\n";
        for (int i = 0; i < argc; ++i)
            cout << argv[i] << ' ';
        cout << endl;
        cout << "</command>\n";
        cout << "-->\n";
        cout << "</html>\n";

        model<doppler_mm_pf> m;
        m.get(cin, cout);

        if (cin.eof() and !cin.fail()) {
            try {
                normal_deviate<mt19937> ndev(seed);
                sym_eig_decomp eig;

                if (m.optimize(parent_number,
                               population_size,
                               step_size,
                               accuracy_goal,
                               stop_generation,
                               trace,
                               ndev, eig, cout))
                    m.put(cout);
            }
            catch (exception &e) {
                cerr << e.what() << endl;

                return 2;
            }
        }
    } else {
        cout << PRG_ID << " " << DOI << endl;
        cout << usemsg << pname << ": " << parmsg << endl;

        return 1;
    }

    return 0;
}
