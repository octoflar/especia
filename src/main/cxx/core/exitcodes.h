/// @file exitcodes.h
/// Application error exit codes.
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
#ifndef ESPECIA_EXITCODES_H
#define ESPECIA_EXITCODES_H

namespace especia {

    /**
     * Application error exit codes.
     */
    class Exit_Codes {
    public:
        /**
         * The optimization stopped due to an underflow of the mutation variance.
         */
        static const int OPTIMIZATION_UNDERFLOW = 1;

        /**
         * The optimization stopped and failed to reach the required accuracy goal.
         */
        static const int OPTIMIZATION_STOPPED = 2;

        /**
         * A logic error occurred.
         */
        static const int LOGIC_ERROR = 10;

        /**
         * A runtime error occurred.
         */
        static const int RUNTIME_ERROR = 11;

        /**
         * An unspecific exception occurred.
         */
        static const int UNSPECIFIC_EXCEPTION = 12;

    private:
        /**
         * Private constructor prevents instantiation.
         */
        Exit_Codes() {

        }
    };

}

#endif // ESPECIA_EXITCODES_H
