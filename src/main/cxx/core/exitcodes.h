/// @file exitcodes.h
/// Application error exit codes.
/// Copyright (c) 2021 Ralf Quast
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
         * The optimization stopped due to an underflow of the mutation variance (exit code = 1).
         */
        static const int optimization_underflow = 00001;

        /**
         * The optimization failed to reach the required accuracy goal within the prescribed
         * number of generations (exit code = 2).
         */
        static const int optimization_stopped = 00002;

        /**
         * A logic error occurred (exit code = 8).
         */
        static const int logic_error = 00010;

        /**
         * A runtime error occurred (exit code = 16).
         */
        static const int runtime_error = 00020;

        /**
         * An unspecific exception occurred (exit code = 64).
         */
        static const int unspecific_exception = 00100;

    private:
        Exit_Codes() = default; // private constructor prevents instantiation
    };

}

#endif // ESPECIA_EXITCODES_H
