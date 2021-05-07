//! @file emod.cxx
//! Program to extract the model definition from an Especia result HTML file

// Copyright (c) 2021. Ralf Quast
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
#include <iostream>

using namespace std;

/**
 * Extracts the model definition from Especia result HTML. Reads from standard
 * input and writes to standard output.
 *
 * @return an exit code.
 *
 * @remark Usage: emod < {result file} [> {target file}]
 */
int main() {
    bool found = false;
    string s;

    while (getline(cin, s)) {
        if (found and not (s == "</model>"
                           // for compatibility with former versions
                           or s == "</job>" or s == "</input>")) {
            cout << s << endl;
        } else {
            found = (s == "<model>"
                     // for compatibility with former versions
                     or s == "<job>" or s == "<input>");
        }
    }

    return 0;
}
