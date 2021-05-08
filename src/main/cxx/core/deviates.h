//! @file deviates.h
//! Function-like class templates for generating various random deviates.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#ifndef ESPECIA_DEVIATES_H
#define ESPECIA_DEVIATES_H

#include <cmath>

#include "base.h"

namespace especia {

   /// A class template to generate random normal deviates.
   ///
   /// The algorithm uses the polar method (e.g. Knuth, 1998,
   /// Sec. 3.4.1, Algorithm P) to generate standard normally distributed random
   /// deviates.
   ///
   /// Further reading:
   ///
   /// D. Knuth (1998).
   ///  *The art of computer programming 2. Seminumerical algorithms.*
   ///  Addison Wesley Longman, ISBN 0-201-89684-2.
   ///
   /// @tparam U The strategy to generate random uniform deviates.
    template<class U>
    class Normal_Deviate {
    public:
        /// Constructs a new instance of this class from a seed.
        ///
        /// @param[in] seed The seed.
        explicit Normal_Deviate(word64 seed = 9600629759793949339ull) : uniform_deviate(seed) {
        }

        /// Constructs a new instance of this class from a uniform deviate.
        ///
        /// @param[in] u The instance of this class to be copied.
        explicit Normal_Deviate(const U &u) : uniform_deviate(u) {
        }

        /// The destructor.
        ~Normal_Deviate() = default;

        /// Returns a normal random number.
        ///
        /// @return a normal random number.
        real operator()() const {
            using std::log;
            using std::sqrt;

            status = !status;

            if (status) {
                real t;

                do {
                    x = 2.0 * uniform_deviate() - 1.0;
                    y = 2.0 * uniform_deviate() - 1.0;
                    t = x * x + y * y;
                } while (t >= 1.0 or t == 0.0);

                t = sqrt(-2.0 * (log(t) / t));
                x *= t;
                y *= t;

                return x;
            } else {
                return y;
            }
        }

    private:
        const U uniform_deviate;

        mutable bool status = false;
        mutable real x = 0.0;
        mutable real y = 0.0;
    };

}

#endif // ESPECIA_DEVIATES_H
