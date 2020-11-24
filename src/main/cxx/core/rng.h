/// @file rng.h
/// Mersenne Twister function-like class template.
/// Copyright (c) 2020 Ralf Quast
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
#ifndef ESPECIA_RNG_H
#define ESPECIA_RNG_H

#include <algorithm>
#include <limits>
#include <valarray>

#include "base.h"

namespace especia {

    /**
     * PCG algorithm  to generate [0,1] uniformly distributed random deviates. Based on
     * Melissa E. O'Neill (2014) and <https://www.pcg-random.org>.
     *
     * Further reading:
     *
     * Melissa E. O'Neill (2014).
     *   *PCG: A Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number Generation.*
     *   <https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf>.
     *
     * @tparam m The multiplier.
     */
    template<word64 m>
    class Pcg {
    public:
        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed The seed.
         * @param[in] selector The sequence selector.
         */
        explicit Pcg(const word64 seed = 9600629759793949339ull, const word64 selector = 7863035247680335341ul) : inc((selector << 1) | 1ull) {
            state = 0ull;
            rand();
            state += seed;
            rand();
        }

        /**
         * The destructor.
         */
        ~Pcg() = default;

        /**
         * Returns a new random number.
         *
         * @return a random number in [0, 1].
         */
        real operator()() const {
            return rand() / real(0xfffffffful);
        }

        word32 rand() const {
            const word64 saved = state;
            state = saved * m + inc;
            const word32 s = (((saved >> 18) ^ saved) >> 27);
            const word32 r = saved >> 59;
            return ((s >> r) | (s << ((-r) & 31u)));
        }

    private:
        /**
         * The increment.
         */
        const word64 inc;
      
        /**
         * The state.
         */
        mutable word64 state;
    };

    /**
     * A predefined PCG algorithm.
     */
    typedef Pcg<6364136223846793005ull> Pcg32;

    /**
     * The Mersenne twister algorithm to generate [0,1] uniformly distributed
     * random deviates.
     *
     * This functor template is based on the 2002/01/26 version coded by Takuji
     * Nishimura and Makoto Matsumoto(Matsumoto and Nishimura, 1998).
     *
     * The notation of template parameters follows Matsumoto and
     * Nishimura (1998, Table 2).
     *
     * Further reading:
     *
     * M. Matsumoto, T. Nishimura (1998).
     *   *Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator.*
     *   ACM Transactions on Modeling and Computer Simulation, 8, 3, ISSN 1049-3301.
     *
     * D. Knuth (1998).
     *   *The art of computer programming 2. Seminumerical algorithms*.
     *   Addison Wesley Longman, ISBN 0-201-89684-2.
     *
     * @tparam w The number of bits in a word.
     * @tparam n The parameter n.
     * @tparam m The parameter m.
     * @tparam r The parameter r.
     * @tparam a The parameter a.
     * @tparam u The parameter u.
     * @tparam s The parameter s.
     * @tparam b The parameter b.
     * @tparam t The parameter t.
     * @tparam c The parameter c.
     * @tparam l The parameter l.
     */
    template<natural w, natural n, natural m, natural r,
            word32 a,
            natural u,
            natural s, word32 b,
            natural t, word32 c,
            natural l>
    class Mersenne_Twister {
    public:
        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed The seed.
         */
        explicit Mersenne_Twister(const word64 seed = 9600629759793949339ull) : words(n) { // NOLINT
            const word32 seeds[] = {word32(seed & 0x00000000ffffffffull), word32(seed & 0xffffffff00000000ull)};
          
            reset(2, seeds);
        }

        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        Mersenne_Twister(const natural seed_count, const word32 seeds[]) : words(n) { // NOLINT
            reset(seed_count, seeds);
        }

        /**
         * The destructor.
         */
        ~Mersenne_Twister() = default;

        /**
         * Returns a new random number.
         *
         * @return a random number in [0, 1].
         */
        real operator()() const {
            using std::numeric_limits;

            // Division by 2^w - 1.
            return rand() * (1.0 / real(numeric_limits<word32>::max() >> (numeric_limits<word32>::digits - w)));
        }

        word32 rand() const {
            if (i == n) {
                for (natural k = 0; k < n - m; ++k) {
                    twist(k + m, k, k + 1);
                }

                for (natural k = n - m; k < n - 1; ++k) {
                    twist(k + m - n, k, k + 1);
                }

                twist(m - 1, n - 1, 0);
                i = 0;
            }

            word32 y = words[i];
            ++i;

            if (u > 0) {
                y ^= (y >> u);
            }

            y ^= (y << s) & b;
            y ^= (y << t) & c;
            y ^= (y >> l);

            return y;
        }
      
    private:
        /**
         * Resets this algorithm.
         *
         * @param[in] seed The seed.
         * @param[in] multiplier A multiplier used by the seeding. Refer to Donald E. Knuth (1997, The Art of Computer
         * Programming, Seminumerical Algorithms, pp. 106) for suitable values.
         */
        void reset(const word32 seed, const word32 multiplier = 1812433253ul) {
            using std::numeric_limits;

            words[0] = seed &
                       (numeric_limits<word32>::max() >> (numeric_limits<word32>::digits - w));
            for (natural k = 1; k < n; ++k) {
                words[k] = (multiplier * (words[k - 1] ^ (words[k - 1] >> (r - 1))) + k) &
                           (numeric_limits<word32>::max() >> (numeric_limits<word32>::digits - w));
            }

            i = n;
        }

        /**
         * Resets this algorithm.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        void reset(const natural seed_count, const word32 seeds[]) {
            using std::max;
            using std::numeric_limits;

            reset(19650218ul);
            i = 1;

            for (natural j = 0, k = max(n, seed_count); k > 0; --k) {
                words[i] = ((words[i] ^ ((words[i - 1] ^ (words[i - 1] >> (r - 1))) * 1664525ul)) + seeds[j] + j) &
                           (numeric_limits<word32>::max() >> (numeric_limits<word32>::digits - w));
                if (++i >= n) {
                    words[0] = words[n - 1];
                    i = 1;
                }
                if (++j >= seed_count) {
                    j = 0;
                }
            }
            for (natural k = n - 1; k > 0; --k) {
                words[i] = ((words[i] ^ ((words[i - 1] ^ (words[i - 1] >> (r - 1))) * 1566083941ul)) - i) &
                           (numeric_limits<word32>::max() >> (numeric_limits<word32>::digits - w));
                if (++i >= n) {
                    words[0] = words[n - 1];
                    i = 1;
                }
            }
            words[0] = (1ul << (w - 1));

            i = n;
        }

        void twist(const natural i, const natural j, const natural k) const {
            using std::numeric_limits;

            words[j] = words[i] ^ (((words[j] & ((numeric_limits<word32>::max()
                    << (numeric_limits<word32>::digits - w + r))
                    >> (numeric_limits<word32>::digits - w))) | (words[k] & (numeric_limits<word32>::max()
                    >> (numeric_limits<word32>::digits - r)))) >> 1);
            if ((words[k] & 1ul) == 1ul) {
                words[j] ^= a;
            }
        }

        mutable std::valarray<word32> words;
        mutable natural i = 0;
    };

    /**
     * A predefined Mersenne twister algorithm.
     */
    typedef Mersenne_Twister<32, 351, 175, 19, 0xe4bd75f5ul, 11, 7, 0x655e5280ul,
            15, 0xffd58000ul, 17> Mt11213a;
    /**
     * A predefined Mersenne twister algorithm.
     */
    typedef Mersenne_Twister<32, 351, 175, 19, 0xccab8ee7ul, 11, 7, 0x31b6ab00ul,
            15, 0xffe50000ul, 17> Mt11213b;
    /**
     * A predefined Mersenne twister algorithm.
     */
    typedef Mersenne_Twister<32, 624, 397, 31, 0x9908b0dful, 11, 7, 0x9d2c5680ul,
            15, 0xefc60000ul, 18> Mt19937;

}

#endif // ESPECIA_RNG_H
