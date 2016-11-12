// Function-like class templates for generating various random deviates
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
