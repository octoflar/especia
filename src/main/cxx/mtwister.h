// Mersenne Twister function-like class template
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
#ifndef RQ_MTWISTER_H
#define RQ_MTWISTER_H

#include <algorithm>
#include <limits>
#include <valarray>

namespace RQ {
    // Mersenne Twister (MT) function-like class template
    template<unsigned w, unsigned n, unsigned m, unsigned r,
            unsigned long a,
            unsigned u,
            unsigned s, unsigned long b,
            unsigned t, unsigned long c,
            unsigned l>
    class mersenne_twister;

    // Provide the MTs tabulated in Matsumoto and Nishimura (1998, Table 2) as
    // predefined types
    typedef mersenne_twister<32, 351, 175, 19, 0xe4bd75f5, 11, 7, 0x655e5280,
            15, 0xffd58000, 17> mt11213a;
    typedef mersenne_twister<32, 351, 175, 19, 0xccab8ee7, 11, 7, 0x31b6ab00,
            15, 0xffe50000, 17> mt11213b;
    typedef mersenne_twister<32, 624, 397, 31, 0x9908b0df, 11, 7, 0x9d2c5680,
            15, 0xefc60000, 18> mt19937;
}

// MT function-like class template for generating [0,1] uniformly distributed
// random deviates, based on the 2002/01/26 version coded by Takuji Nishimura
// and Makoto Matsumoto
//
// The notation of the template parameters follows Matsumoto
// and Nishimura (1998,Table 2).
template<unsigned w, unsigned n, unsigned m, unsigned r,
        unsigned long a,
        unsigned u,
        unsigned s, unsigned long b,
        unsigned t, unsigned long c,
        unsigned l>
class RQ::mersenne_twister {
public:
    // Initialize the MT using an LCG with modulus 2^w. Knuth (1998,
    // Sec. 3.3.4, Table 1) gives some reasonably good ones.
    mersenne_twister(unsigned long seed = 5489, unsigned long multiplier = 1812433253)
            : x(n) {
        reset(seed, multiplier);
    }

    // Initialize the MT using an array of seed values
    mersenne_twister(const unsigned long key[], unsigned key_size)
            : x(n) {
        reset(key, key_size);
    }

    ~mersenne_twister() {
    };

    double operator()() {
        using std::numeric_limits;

        return rand() * (1.0 / double(numeric_limits<word>::max() >> (numeric_limits<word>::digits - w)));
            // division by 2^w - 1
    }

    void reset(unsigned long seed = 5489, unsigned long multiplier = 1812433253) {
        using std::max;
        using std::numeric_limits;

        x[0] = seed & (numeric_limits<word>::max() >> (numeric_limits<word>::digits - w));
        for (unsigned k = 1; k < n; ++k) {
            x[k] = (multiplier * (x[k - 1] ^ (x[k - 1] >> (r - 1))) + k) &
                   (numeric_limits<word>::max() >> (numeric_limits<word>::digits - w));
        }

        i = n;
    }

    void reset(const unsigned long key[], unsigned key_size) {
        using std::max;
        using std::numeric_limits;

        reset(19650218);
        i = 1;

        for (unsigned j = 0, k = max(n, key_size); k > 0; --k) {
            x[i] = ((x[i] ^ ((x[i - 1] ^ (x[i - 1] >> (r - 1))) * 1664525)) + key[j] + j) &
                   (numeric_limits<word>::max() >> (numeric_limits<word>::digits - w));
            if (++i >= n) {
                x[0] = x[n - 1];
                i = 1;
            }
            if (++j >= key_size)
                j = 0;
        }
        for (unsigned k = n - 1; k > 0; --k) {
            x[i] = ((x[i] ^ ((x[i - 1] ^ (x[i - 1] >> (r - 1))) * 1566083941)) - i) &
                   (numeric_limits<word>::max() >> (numeric_limits<word>::digits - w));
            if (++i >= n) {
                x[0] = x[n - 1];
                i = 1;
            }
        }
        x[0] = (1ul << (w - 1));
        // MSB is 1, assuring a non-zero initial array

        i = n;
    }

private:
    typedef unsigned long word;

    unsigned long rand() {
        if (i == n) {
            for (unsigned k = 0; k < n - m; ++k)
                twist(k + m, k, k + 1);

            for (unsigned k = n - m; k < n - 1; ++k)
                twist(k + m - n, k, k + 1);

            twist(m - 1, n - 1, 0);
            i = 0;
        }

        unsigned long y = x[i];
        ++i;

        if (u > 0)
            y ^= (y >> u);

        y ^= (y << s) & b;
        y ^= (y << t) & c;
        y ^= (y >> l);

        return y;
    }

    void twist(unsigned i, unsigned j, unsigned k) {
        using std::numeric_limits;

        x[j] = x[i] ^ (((x[j] & ((numeric_limits<word>::max() << (numeric_limits<word>::digits - w + r))
                >> (numeric_limits<word>::digits - w))) |
                        (x[k] & (numeric_limits<word>::max() >> (numeric_limits<word>::digits - r)))) >> 1);
        if ((x[k] & 1ul) == 1ul)
            x[j] ^= a;
    }

    std::valarray<unsigned long> x;
    unsigned i;
};

#endif // RQ_MTWISTER_H

// References
//
// D. Knuth (1998)
//   The art of computer programming 2. Seminumerical algorithms
//   Addison Wesley Longman, ISBN 0-201-89684-2
//
// M. Matsumoto, T. Nishimura (1998)
//   Mersenne Twister: A 623-dimensionally equidistributed uniform
//   pseudorandom number generator
//   ACM Transactions on Modeling and Computer Simulation, 8, 3,
//   ISSN 1049-3301
//
