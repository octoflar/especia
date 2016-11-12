// Function-like class templates for generating various random deviates
// Copyright (c) 2016, Ralf Quast
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef RQ_RANDEV_H
#define RQ_RANDEV_H

#include <cmath>

namespace RQ {
    // Function-like class template for generating normal deviates
    template<class uniform_deviate>
    class normal_deviate;
}

// Function-like class template using the polar method (e.g. Knuth, 1998,
// Sec. 3.4.1, Algorithm P) to generate standard normally distributed random
// deviates
template<class uniform_deviate>
class RQ::normal_deviate {
public:
    normal_deviate(unsigned long seed = 5489);
    normal_deviate(uniform_deviate& udev);
    ~normal_deviate();

    double operator()();
    void reset(unsigned long seed = 5489);
    void reset(uniform_deviate& udev);

private:
    uniform_deviate udev;
    bool gen_xy;
    double x, y;
};

template<class uniform_deviate>
RQ::normal_deviate<uniform_deviate>::normal_deviate(unsigned long seed)
    :   udev(seed), gen_xy(false)
{
}

template<class uniform_deviate>
RQ::normal_deviate<uniform_deviate>::normal_deviate(uniform_deviate& u)
    :   udev(u), gen_xy(false)
{
}

template<class uniform_deviate>
RQ::normal_deviate<uniform_deviate>::~normal_deviate()
{
}

template<class uniform_deviate>
double
RQ::normal_deviate<uniform_deviate>::operator()()
{
    using std::log;
    using std::sqrt;

    gen_xy = !gen_xy;

    if (gen_xy) {
        double t;

        do {
            x = 2.0 * udev() - 1.0;
            y = 2.0 * udev() - 1.0;
            t = x * x + y * y;
        } while (t >= 1.0 or t == 0.0);

        t = sqrt(-2.0 * (log(t) / t));
        x *= t;
        y *= t;

        return x;
    } else
        return y;
}

template<class uniform_deviate>
void
RQ::normal_deviate<uniform_deviate>::reset(unsigned long seed)
{
    udev.reset(seed);
    gen_xy = false;
}

template<class uniform_deviate>
void
RQ::normal_deviate<uniform_deviate>::reset(uniform_deviate& u)
{
    udev = u;
    gen_xy = false;
}

#endif // RQ_RANDEV_H

// References
//
// D. Knuth (1998)
//   The art of computer programming 2. Seminumerical algorithms
//   Addison Wesley Longman, ISBN 0-201-89684-2
