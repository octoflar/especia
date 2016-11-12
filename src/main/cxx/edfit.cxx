// Darwin (Doppler profile version)
// Copyright (c) 2016, Ralf Quast
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <exception>
#include <fstream>
#include "mtwister.h"
#include "randev.h"
#include "model.h"
#include "symeig.h"

using namespace RQ;
using namespace std;

const char usemsg[] = "Usage: ";
const char parmsg[] = "SEED PARENTS POPULATION INISTEP ACCURACY STOPGEN TRACE";

int main(int argc, char* argv[])
{
    const char* pname = argv[0];

    if (argc == 8) {
        const long   seed            = atol(argv[1]);
        const int    parent_number   = atol(argv[2]);
        const int    population_size = atoi(argv[3]);
        const double step_size       = atof(argv[4]);
        const double accuracy_goal   = atof(argv[5]);
        const long   stop_generation = atol(argv[6]);
        const int    trace           = atoi(argv[7]);

        model<gauss_pf> m;
        m.get(cin);

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
            catch (exception& e) {
                cerr << e.what() << endl;

                return 0;
            }
        }
    } else {
        cerr << usemsg << pname << ": " << parmsg << endl;

        return 1;
    }

    return 0;
}
