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
#include <sstream>
#include <string>


/**
 * The base class to be inherited by all unit-level tests.
 */
class Unit_Test {
public:
    /**
     * The type of exception thrown when an assertion fails.
     */
    class Assertion_Error : public std::exception {
    public:
        /**
         * The constructor.
         *
         * @param what A description of the failed assertion.
         */
        Assertion_Error(std::string what) : exception(), what_happened(what) {
        }

        /**
         * The destructor.
         */
        virtual ~Assertion_Error() {
        }

        /**
         * Returns a description of the failed assertion.
         *
         * @return the description of the failed assertion.
         */
        virtual const char* what() const noexcept {
            return what_happened.c_str();
        }

    private:
        /**
         * The description of the failed assertion.
         */
        std::string what_happened;
    };

    /**
     * The constructor.
     */
    Unit_Test() {

    }

    /**
     * The destructor.
     */
    virtual ~Unit_Test() {

    }

    /**
     * Runs the testsuite.
     *
     * @throw an @c Assertion_Error when an assertion fails.
     */
    void run_testsuite() throw(Assertion_Error) {
        before_all();
        run_all();
        after_all();
    }

protected:
    /**
     * Method called before any test case will be executed.
     */
    virtual void before_all() {

    }

    /**
     * Method called after all test cases have been executed.
     */
    virtual void after_all() {

    }

    /**
     * Method called before each test case.
     */
    virtual void before() {

    }

    /**
     * Method called after each test case.
     */
    virtual void after() {

    }

    /**
     * Runs all test cases.
     */
    virtual void run_all() = 0;

    /**
     * Runs a test case.
     *
     * @tparam C The test class.
     * @tparam T The test case type.
     *
     * @param c A pointer to the test class.
     * @param t A pointer to the test case.
     */
    template<class C, class T>
    void run(C c, T t) {
        before();
        (c->*t)();
        after();
    }

    /**
     * Asserts equality of two values.
     *
     * @tparam T The value type.
     *
     * @param name The assertion name.
     * @param expected The expected value.
     * @param actual The actual value.
     */
    template<class T>
    void assert_equals(const T &expected, const T &actual,
                       const std::string &name = "unnamed assertion") const throw(Assertion_Error) {
        using std::cerr;
        using std::endl;
        using std::stringstream;

        if (actual == expected) { // NaN safe
            conditional_success_message(name);
        } else {
            stringstream what;
            what << "Failed: " << name << "\n";
            what << "    Expected result: " << expected << "\n";
            what << "    Actual   result: " << actual;
            if (message_on_failure) {
                cerr << what.str() << endl;
            }
            conditional_throw(what.str());
        }
    }

    /**
     * Asserts equality of two values with some (absolute) tolerance.
     *
     * @tparam T The value type.
     *
     * @param expected The expected value.
     * @param actual The actual value.
     * @param tolerance The (absolute) tolerance.
     * @param name The assertion name.
     */
    template<class T>
    void assert_equals(const T &expected, const T &actual, const T &tolerance,
                       const std::string &name = "unnamed assertion") const throw(Assertion_Error) {
        using std::abs;
        using std::cerr;
        using std::endl;
        using std::stringstream;

        if (abs(actual - expected) < tolerance) { // NaN safe
            conditional_success_message(name);
        } else {
            stringstream what;
            what << "Failed: " << name << "\n";
            what << "    Expected result: " << expected << "\n";
            what << "    Actual   result: " << actual << "\n";
            what << "    Expected tolerance: " << tolerance;
            if (message_on_failure) {
                cerr << what.str() << endl;
            }
            conditional_throw(what.str());
        }
    }

private:
    /**
     * Throws an assertion error when an assertion is failed, if requested.
     *
     * @param what A description of the failed assertion.
     * @throw an @c Assertion_Error if requested.
     */
    void conditional_throw(const std::string &what) const throw(Assertion_Error) {
        if (throw_on_failure) {
            throw Assertion_Error(what);
        }
    }

    /**
     * Issues a success message when an assertion is passed, if requested.
     *
     * @param name The name of the passed assertion.
     */
    void conditional_success_message(const std::string &name) const {
        using std::cout;
        using std::endl;

        if (message_on_success) {
            cout << "Passed: " << name << endl;
        }
    }

    /**
     * Issue a message when an assertion is failed?
     */
    bool message_on_failure = false;

    /**
     * Issue a message when an assertion is passed?
     */
    bool message_on_success = true;

    /**
     * Throw an exception when an assertion is failed?
     */
    bool throw_on_failure = true;
};

#endif // ESPECIA_UNITTEST_H_H
