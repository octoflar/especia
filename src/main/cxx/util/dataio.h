// Data input and output procedures
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

#ifndef ESPECIA_DATAIO_H
#define ESPECIA_DATAIO_H

#include <iostream>
#include <valarray>

namespace especia {

    /**
     * Reads spectroscopic data from an input stream (two-column format).
     *
     * @param is[in,out] The input stream.
     * @param x[out] The wavelength data (previous content is overwritten).
     * @param y[out] The spectral flux or uncertainty data (previous content is overwritten).
     * @param skip[in] The number of leading lines to skip (optional).
     *
     * @return the input stream.
     */
    std::istream &
    get(std::istream &is, std::valarray<double> &x, std::valarray<double> &y, int skip = 0);

    /**
     * Reads spectroscopic data from an input stream.
     *
     * @param is[in,out] The input stream.
     * @param x[out] The wavelength data (previous content is overwritten).
     * @param y[out] The spectral flux data (previous content is overwritten).
     * @param z[out] The spectral flux uncertainty data (previous content is overwritten).
     * @param skip[in] The number of leading lines to skip (optional).
     *
     * @return the input stream.
     */
    std::istream &
    get(std::istream &is, std::valarray<double> &x, std::valarray<double> &y, std::valarray<double> &z, int skip = 0);

    /**
     * Writes spectroscopic data to an output stream.
     *
     * @param os[in,out] The output stream.
     * @param x[in] The wavelength data.
     * @param y[in] The spectral flux data.
     * @param z[in] The spectral flux uncertainty data.
     *
     * @return the output stream.
     */
    std::ostream &
    put(std::ostream &os, const std::valarray<double> &x, const std::valarray<double> &y,
        const std::valarray<double> &z);

}

#endif // ESPECIA_DATAIO_H
