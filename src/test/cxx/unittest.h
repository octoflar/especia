/// @file unittest.h
/// Simple unit testing framework
/// Copyright (c) 2017 Ralf Quast
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
#ifndef ESPECIA_UNITTEST_H_H
#define ESPECIA_UNITTEST_H_H

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using std::runtime_error;
using std::string;


class Assertions {
public:

    Assertions() {

    }

    ~Assertions() {

    }

    template<class T>
    void assert_equals(const string &name, const T &expected, const T &actual, const T &tolerance) const {
        using std::abs;
        using std::cerr;
        using std::endl;

        if (abs(actual - expected) > tolerance) {
            if (message_on_failure) {
                cerr << "Assertions: assertion '" << name << "' failed" << std::endl;
                cerr << "    Expected result: " << expected << endl;
                cerr << "    Actual   result: " << actual << endl;
                cerr << "    Expected tolerance: " << tolerance << endl;
            }
            conditional_stop(name);
        } else {
            conditional_success_message(name);
        }
    }

private:
    void conditional_stop(const string &name) const throw(runtime_error) {
        if (stop_on_failure) {
            throw runtime_error("Assertions: assertion '" + name + "' failed");
        }
    }

    void conditional_success_message(const string &name) const {
        using std::cout;
        using std::endl;

        if (message_on_success) {
            cout << "Assertions: assertion '" << name << "' passed" << endl;
        }
    }

    bool message_on_failure = true;
    bool message_on_success = true;
    bool stop_on_failure = true;
};

#endif // ESPECIA_UNITTEST_H_H
