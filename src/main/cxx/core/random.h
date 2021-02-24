/// @file random.h
/// Function-like class templates to generate random numbers.
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
#ifndef ESPECIA_RANDOM_H
#define ESPECIA_RANDOM_H

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
     * @tparam mult The multiplier.
     */
    template<word64 mult>
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
         * Returns a new real-valued  random number in the interval  [0, 1].
         *
         * @return a real-valued random number in [0, 1].
         */
        real operator()() const {
            return rand() / real(0xfffffffful);
        }

        /**
         * Returns a new random word.
         *
         * @return a random word.
         */
        word32 rand() const {
            const word64 saved = state;
            state = saved * mult + inc;
            const word32 s = (((saved >> 18) ^ saved) >> 27);
            const word32 r = saved >> 59;
            return ((s >> r) | (s << ((-r) & 31ul)));
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
     * The PCG-XSH-RR with 64-bit state and 32-bit output.
     */
    typedef Pcg<6364136223846793005ull> Pcg_32;


    /**
     * A maximally equidistributed F2-linear generator (MELG). This MELG is formally
     * designated 'MELG19937-64'. It has a state of 2,496 bytes and yields a 64-bit
     * output word.
     *
     * Further reading:
     *
     * S. Harase and T. Kimoto (2018).
     * Implementing 64-bit maximally equidistributed F2-linear generators with Mersenne prime period.
     * ACM Transactions on Mathematical Software, 44, 3, 30.
     * <http://doi.acm.org/10.1145/3159444>, <http://arxiv.org/abs/1505.06582>
     *
     * @tparam w The the number of bits in a word.
     */
    template<natural w>
    class Melg {
    public:
        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed The seed.
         */
        explicit Melg(const word64 seed) : state(N + 1) { // NOLINT
            reset(seed);
        }

        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        Melg(const natural seed_count, const word64 seeds[]) : state(N + 1) { // NOLINT
            reset(seed_count, seeds);
        }

        /**
         * The destructor.
         */
        ~Melg() = default;

        /**
         * Returns a new real-valued  random number in the interval  [0, 1].
         *
         * @return a real-valued random number in [0, 1].
         */
        real operator()() const {
            using std::numeric_limits;
            // the effective maximum mantissa value for a real number
            const real max_real_mantissa =
                    numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - (w < numeric_limits<real>::digits ? w : numeric_limits<real>::digits));

            return (w < numeric_limits<real>::digits ? rand() : rand() >> (w - numeric_limits<real>::digits)) * (1.0 / max_real_mantissa);
        }

        /**
         * Returns a new  random word.
         *
         * @return a random word.
         */
        word64 rand() const {
            word64 next = 0ull;

            switch (state_c) {
                case 1:
                    next = twist1(state_i, state_i + 1);
                    twist2(next, state_i + M);
                    next = twist3(next, state_i, state_i + L);
                    state_i++;
                    if (state_i == N - M) {
                        state_c = 2;
                    }
                    break;
                case 2:
                    next = twist1(state_i, state_i + 1);
                    twist2(next, state_i + M - N);
                    next = twist3(next, state_i, state_i + L);
                    state_i++;
                    if (state_i == N - L) {
                        state_c = 3;
                    }
                    break;
                case 3:
                    next = twist1(state_i, state_i + 1);
                    twist2(next, state_i + M - N);
                    next = twist3(next, state_i, state_i - (N - L));
                    state_i++;
                    if (state_i == N - 1) {
                        state_c = 4;
                    }
                    break;
                case 4:
                    next = twist1(N - 1, 0);
                    twist2(next, M - 1);
                    next = twist3(next, N - 1, state_i - (N - L));
                    state_i = 0;
                    state_c = 1;
                    break;
            }

            return next;
        }

    private:
        /**
         * Resets this algorithm.
         *
         * @param[in] seed The seed.
         */
        void reset(const word64 seed) {
            state[0] = seed;

            for (state_i = 1; state_i < N + 1; state_i++) {
                state[state_i] = (state[state_i - 1] ^ (state[state_i - 1] >> 62)) * 6364136223846793005ull + state_i;
            }

            state_i = 0;
            state_c = 1;
        }

        /**
         * Resets this algorithm.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        void reset(const natural seed_count, const word64 seeds[]) {
            using std::max;

            reset(19650218ull);

            natural i = 1;
            natural j = 0;
            for (natural k = max(N, seed_count); k > 0; k--) {
                state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 62)) * 3935559000370003845ull)) + seeds[j] + j;
                i++;
                j++;
                if (i >= N) {
                    state[0] = state[N - 1];
                    i = 1;
                }
                if (j >= seed_count) {
                    j = 0;
                }
            }
            for (natural k = N - 1; k > 0; k--) {
                state[i] = (state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 62)) * 2862933555777941757ull)) - i;
                i++;
                if (i >= N) {
                    state[0] = state[N - 1];
                    i = 1;
                }
            }
            state[N] = (state[N] ^ ((state[N - 1] ^ (state[N - 1] >> 62)) * 2862933555777941757ull)) - N;
            state[0] = (state[0] | (1ull << 63));
            state_i = 0;
            state_c = 1;
        }

        word64 twist1(const natural i1, const natural i2) const {
            return (state[i1] & 0xffffffff80000000ull) | (state[i2] & 0x7fffffffull);
        }

        void twist2(const word64 l1, const natural i1) const {
            state[N] = (l1 >> 1) ^ (((l1 & 1ull) != 0ull) ? 0x5c32e06df730fc42ull : 0ull) ^ state[i1] ^ (state[N] ^ (state[N] << 23));
        }

        word64 twist3(const word64 l1, const natural i1, const natural i2) const {
            state[i1] = l1 ^ (state[N] ^ (state[N] >> 33));
            return state[i1] ^ (state[i1] << 16) ^ (state[i2] & 0x6aede6fd97b338ecull);
        }

        mutable std::valarray<word64> state;
        mutable natural state_i;
        mutable natural state_c;

        static const natural L = 19;
        static const natural M = 81;
        static const natural N = 311;
    };

    /**
     * The MELG19937-64 with 2,496 bytes of state and 64-bit output.
     */
    typedef Melg<64> Melg19937_64;


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
     * @tparam d The parameter d.
     * @tparam s The parameter s.
     * @tparam b The parameter b.
     * @tparam t The parameter t.
     * @tparam c The parameter c.
     * @tparam l The parameter l.
     * @tparam mult1 The first multiplier (used for initialisation).
     * @tparam mult2 The second multiplier (used for initialisation).
     * @tparam mult3 The third multiplier (used for initialisation).
     */
    template<natural w, natural n, natural m, natural r,
            word64 a,
            natural u, word64 d,
            natural s, word64 b,
            natural t, word64 c,
            natural l,
            word64 mult1,
            word64 mult2,
            word64 mult3>
    class Mersenne_Twister {
    public:
        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed The seed.
         */
        explicit Mersenne_Twister(const word64 seed = 9600629759793949339ull) : words(n) { // NOLINT
            const word64 seeds[] = {seed & 0x00000000ffffffffull, seed & 0xffffffff00000000ull};

            reset(2, seeds);
        }

        /**
         * Constructs a new instance of this functor.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        Mersenne_Twister(const natural seed_count, const word64 seeds[]) : words(n) { // NOLINT
            reset(seed_count, seeds);
        }

        /**
         * The destructor.
         */
        ~Mersenne_Twister() = default;

        /**
         * Returns a new real-valued  random number in the interval  [0, 1].
         *
         * @return a real-valued random number in [0, 1].
         */
        real operator()() const {
            using std::numeric_limits;
            // the effective maximum mantissa value for a real number
            const real max_real_mantissa =
                numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - (w < numeric_limits<real>::digits ? w : numeric_limits<real>::digits));
            
            return (w < numeric_limits<real>::digits ? rand() : rand() >> (w - numeric_limits<real>::digits)) * (1.0 / max_real_mantissa);
        }

        /**
         * Returns a new  random word.
         *
         * @return a random word.
         */
        word64 rand() const {
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

            word64 y = words[i];
            ++i;

            y ^= (y >> u) & d;
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
         */
        void reset(const word64 seed) {
            using std::numeric_limits;

            words[0] = seed &
                       (numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - w));
            for (natural k = 1; k < n; ++k) {
                words[k] = (mult1 * (words[k - 1] ^ (words[k - 1] >> (w - 2))) + k) &
                           (numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - w));
            }

            i = n;
        }

        /**
         * Resets this algorithm.
         *
         * @param[in] seed_count The number of seeds.
         * @param[in] seeds The seeds.
         */
        void reset(const natural seed_count, const word64 seeds[]) {
            using std::max;
            using std::numeric_limits;

            reset(19650218ull);
            i = 1;

            for (natural j = 0, k = max(n, seed_count); k > 0; --k) {
                words[i] = ((words[i] ^ ((words[i - 1] ^ (words[i - 1] >> (w - 2))) * mult2)) + seeds[j] + j) &
                           (numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - w));
                if (++i >= n) {
                    words[0] = words[n - 1];
                    i = 1;
                }
                if (++j >= seed_count) {
                    j = 0;
                }
            }
            for (natural k = n - 1; k > 0; --k) {
                words[i] = ((words[i] ^ ((words[i - 1] ^ (words[i - 1] >> (w - 2))) * mult3)) - i) &
                           (numeric_limits<word64>::max() >> (numeric_limits<word64>::digits - w));
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

            words[j] = words[i] ^ (((words[j] & ((numeric_limits<word64>::max()
                    << (numeric_limits<word64>::digits - w + r))
                    >> (numeric_limits<word64>::digits - w))) | (words[k] & (numeric_limits<word64>::max()
                    >> (numeric_limits<word64>::digits - r)))) >> 1);
            if ((words[k] & 1ull) == 1ull) {
                words[j] ^= a;
            }
        }

        mutable std::valarray<word64> words;
        mutable natural i = 0;
    };

    /**
     * The MT-1121A-32.
     */
    typedef Mersenne_Twister<32, 351, 175, 19, 0xe4bd75f5ull, 11, 0xffffffffull, 7, 0x655e5280ull,
            15, 0xffd58000ull, 17, 1812433253ull, 1664525ull, 1566083941ull> Mt11213a_32;
    /**
     * The MT-1121B-32.
     */
    typedef Mersenne_Twister<32, 351, 175, 19, 0xccab8ee7ull, 11, 0xffffffffull, 7, 0x31b6ab00ull,
            15, 0xffe50000ull, 17, 1812433253ull, 1664525ull, 1566083941ull> Mt11213b_32;
    /**
     * The MT-19937-32.
     */
    typedef Mersenne_Twister<32, 624, 397, 31, 0x9908b0dfull, 11, 0xffffffffull, 7, 0x9d2c5680ull,
            15, 0xefc60000ull, 18, 1812433253ul, 1664525ull, 1566083941ull> Mt19937_32;

    /**
     * The MT-19937-64.
     */
    typedef Mersenne_Twister<64, 312, 156, 31, 0xB5026F5AA96619E9ull, 29, 0x5555555555555555ull, 17, 0x71D67FFFEDA60000ull,
            37, 0xFFF7EEE000000000ull, 43, 6364136223846793005ull, 3935559000370003845ull, 2862933555777941757ull> Mt19937_64;
}

#endif // ESPECIA_RANDOM_H
