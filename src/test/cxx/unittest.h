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
#include <exception>
#include <iostream>
#include <string>


class Unit_Test {
public:
    class Assertion_Error : public std::exception {
    public:
        Assertion_Error(std::string description) : exception(), text(description) {
        }

        virtual ~Assertion_Error() {
        }

        virtual const char* what() const _NOEXCEPT {
            return text.c_str();
        }

    private:
        std::string text;
    };

    virtual ~Unit_Test() {

    }

    void run_testsuite() throw(Assertion_Error) {
        before_all();
        run_all();
        after_all();
    }

protected:
    Unit_Test() {

    }

    virtual void before_all() {

    }

    virtual void after_all() {

    }

    virtual void before() {

    }

    virtual void after() {

    }

    template<class C, class T>
    void run(C class_pointer, T test_pointer) {
        before();
        (class_pointer->*test_pointer)();
        after();
    }

    virtual void run_all() = 0;

    template<class T>
    void assert_equals(const std::string &name,
                       const T &expected, const T &actual) const throw(Assertion_Error) {
        using std::cerr;
        using std::endl;

        if (actual == expected) { // NaN safe
            conditional_success_message(name);
        } else {
            if (message_on_failure) {
                cerr << "Failed: " << name << "\n";
                cerr << "    Expected result: " << expected << "\n";
                cerr << "    Actual   result: " << actual << endl;
            }
            conditional_stop(name);
        }
    }

    template<class T>
    void assert_equals(const std::string &name,
                       const T &expected, const T &actual, const T &tolerance) const throw(Assertion_Error) {
        using std::abs;
        using std::cerr;
        using std::endl;

        if (abs(actual - expected) < tolerance) { // NaN safe
            conditional_success_message(name);
        } else {
            if (message_on_failure) {
                cerr << "Failed: " << name << "\n";
                cerr << "    Expected result: " << expected << "\n";
                cerr << "    Actual   result: " << actual << "\n";
                cerr << "    Expected tolerance: " << tolerance << endl;
            }
            conditional_stop(name);
        }
    }

private:
    void conditional_stop(const std::string &name) const throw(Assertion_Error) {
        if (stop_on_failure) {
            throw Assertion_Error("Failed: " + name);
        }
    }

    void conditional_success_message(const std::string &name) const {
        using std::cout;
        using std::endl;

        if (message_on_success) {
            cout << "Passed: " << name << endl;
        }
    }

    bool message_on_failure = true;
    bool message_on_success = true;
    bool stop_on_failure = true;
};

#endif // ESPECIA_UNITTEST_H_H
