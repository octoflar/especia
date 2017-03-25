/// @file integrator.h
/// Numerical integration.
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
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "base.h"

namespace especia {

    /**
     * Numerical integration by means of recursive monotone stable quadrature
     * formulas.
     *
     * Further reading:
     *
     * Favati, P.; Lotti, G.; and Romani, F. (1991).
     *   *Interpolary Integration Formulas for Optimal Composition.*
     *   ACM Trans. Math. Software 17, 207-217.
     *   https://doi.org/10.1145/108556.108571
     * Favati, P.; Lotti, G.; and Romani, F. (1991).
     *   *Algorithm 691: Improving QUADPACK Automatic Integration Routines.*
     *   ACM Trans. Math. Software 17, 218-232.
     *   https://doi.org/10.1145/108556.108580
     *
     * @tparam T The number type used for quadrature calculations. The quadrature weight
     * and abscissa values are defined with a precision of 48 decimal digits.
     *
     * @todo - work in progress
     */
    template<class T = R_type>
    class Integrator {
    public:

        /**
         * The recursive monotone stable quadrature formulas.
         */
        enum Formula {
            /**
             * The formula for integration with 13 quadrature points.
             */
                    Q13,
            /**
             * The formula for integration with 19 quadrature points.
             */
                    Q19,
            /**
             * The formula for integration with 27 quadrature points.
             */
                    Q27,
            /**
             * The formula for integration with 41 quadrature points.
             */
                    Q41
        };

        /**
         * Constructs a new integrator, based on the formula supplied as argument.
         *
         * @param q
         */
        Integrator(Formula q = Q27) : q(q) {
        }

        /**
         * Destructor.
         */
        ~Integrator() {

        }

        /**
         * Computes the value of the integral of a function, with the limits supplied as argument.
         *
         * @tparam F The function type.
         *
         * @param[in] f The function.
         * @param[in] a The lower limit of integration.
         * @param[in] b The upper limit of integration.
         * @return the value of the integral.
         */
        template<class F>
        T integrate(const F &f, T a, T b) {
            const size_t m = mw[q];
            const size_t n = nw[q];

            const T c = T(0.5) * (a + b);
            const T h = T(0.5) * (b - a);

            T result = T(0.0);

            for (size_t i = 0; i < n; ++i) {
                result += (f(c - h * xi[i]) + f(c + h * xi[i])) * wi[m + i];
            }

            return result * h;
        }

    private:
        friend class Interval;

        class Interval {
        public:

            Interval(T a_in, T b_in)
                    : a(a_in), b(b_in), c(T(0.5) * (a + b)), h(T(0.5) * (b - a)), ym(21), yn(21) {
            }

            ~Interval() {

            }

            template<class F>
            Interval lower_half(const F &f) const {
                Interval interval(a, c);

                interval.yn[0] = ym[2];
                interval.yn[1] = ym[7];
                interval.yn[2] = ym[1];
                interval.yn[4] = f(interval.c + xi[4] * interval.h);
                interval.yn[5] = f(interval.c + xi[5] * interval.h);
                interval.yn[6] = ym[0];
                interval.ym[0] = ym[2];
                interval.ym[1] = ym[8];
                interval.ym[2] = ym[3];
                interval.ym[3] = ym[4];
                interval.ym[4] = ym[5];
                interval.ym[5] = ym[9];
                interval.ym[6] = ym[6];
                if (m > 10) {
                    interval.yn[3] = ym[10];
                    interval.ym[7] = ym[11];
                    interval.ym[8] = ym[12];
                    interval.ym[9] = ym[13];
                    if (m > 14) {
                        interval.yn[7] = ym[15];
                        interval.yn[8] = ym[14];
                        interval.yn[9] = f(interval.c + xi[9] * interval.h);
                        interval.yn[10] = ym[16];
                        interval.ym[10] = ym[17];
                        interval.ym[11] = ym[18];
                        interval.ym[12] = ym[19];
                        interval.ym[13] = ym[20];
                        interval.n = 11;
                        interval.m = 14;
                    } else {
                        interval.n = 7;
                        interval.m = 10;
                    }
                } else {
                    interval.yn[3] = f(interval.c + xi[3] * interval.h);
                    interval.n = 7;
                    interval.m = 7;
                }

                return interval;
            }

            template<class F>
            Interval upper_half(const F &f) const {
                Interval interval(c, b);

                interval.ym[0] = yn[2];
                interval.ym[1] = yn[7];
                interval.ym[2] = yn[1];
                interval.ym[4] = f(interval.c - xi[4] * interval.h);
                interval.ym[5] = f(interval.c - xi[5] * interval.h);
                interval.ym[6] = yn[0];
                interval.yn[0] = yn[2];
                interval.yn[1] = yn[8];
                interval.yn[2] = yn[3];
                interval.yn[3] = yn[4];
                interval.yn[4] = yn[5];
                interval.yn[5] = yn[9];
                interval.yn[6] = yn[6];
                if (n > 10) {
                    interval.ym[3] = yn[10];
                    interval.yn[7] = yn[11];
                    interval.yn[8] = yn[12];
                    interval.yn[9] = yn[13];
                    if (n > 14) {
                        interval.ym[7] = yn[15];
                        interval.ym[8] = yn[14];
                        interval.ym[9] = f(interval.c - xi[9] * interval.h);
                        interval.ym[10] = yn[16];
                        interval.yn[10] = yn[17];
                        interval.yn[11] = yn[18];
                        interval.yn[12] = yn[19];
                        interval.yn[13] = yn[20];
                        interval.m = 11;
                        interval.n = 14;
                    } else {
                        interval.m = 7;
                        interval.n = 10;
                    }
                } else {
                    interval.ym[3] = f(interval.c - xi[3] * interval.h);
                    interval.m = 7;
                    interval.n = 7;
                }

                return interval;
            }

            const T a;
            const T b;
            const T c;
            const T h;

            std::valarray<T> ym;
            std::valarray<T> yn;

            size_t m = 0;
            size_t n = 0;

            T result = T(0.0);
            T abserr = T(0.0);
        };

        /**
         * The selected quadrature formula.
         */
        const Formula q;

        /**
         * The quadrature abscissa values.
         */
        static const T xi[];

        /**
         * The quadrature weights.
         */
        static const T wi[];

        /**
         * The start indices into the quadrature weights.
         */
        static const size_t mw[];

        /**
         * The number of quadrature weights.
         */
        static const size_t nw[];
    };

    template<class T>
    const T Integrator<T>::xi[] = {
            // abscissas for Q13
            0.0000000,
            0.2500000,
            0.5000000,
            0.7500000,
            0.8750000,
            0.9375000,
            1.0000000,
            // additional abscissas for Q19, Q27 and Q41
            0.3750000,
            0.6250000,
            0.9687500,
            // additional abscissas for Q27 and Q41
            0.1250000,
            0.6875000,
            0.8125000,
            0.9843750,
            // additional abscissas for Q41
            0.1875000,
            0.3125000,
            0.4375000,
            0.5625000,
            0.8437500,
            0.9062500,
            0.9921875
    };

    template<class T>
    const T Integrator<T>::wi[] = {
            // weights for Q13
            1.303262173284849021810473057638590518409112513421E-01,
            2.390632866847646220320329836544615917290026806242E-01,
            2.630626354774670227333506083741355715758124943143E-01,
            2.186819313830574175167853094864355208948886875898E-01,
            2.757897646642836865859601197607471574336674206700E-02,
            1.055750100538458443365034879086669791305550493830E-01,
            1.571194260595182254168429283636656908546309467968E-02,
            // weights for Q19
            1.298751627936015783241173611320651866834051160074E-01,
            2.249996826462523640447834514709508786970828213187E-01,
            1.680415725925575286319046726692683040162290325505E-01,
            1.415567675701225879892811622832845252125600939627E-01,
            1.006482260551160175038684459742336605269707889822E-01,
            2.510604860724282479058338820428989444699235030871E-02,
            9.402964360009747110031098328922608224934320397592E-03,
            5.542699233295875168406783695143646338274805359780E-02,
            9.986735247403367525720377847755415293097913496236E-02,
            4.507523056810492466415880450799432587809828791196E-02,
            // weights for Q27
            6.300942249647773931746170540321811473310938661469E-02,
            1.261383225537664703012999637242003647020326905948E-01,
            1.273864433581028272878709981850307363453523117880E-01,
            8.576500414311820514214087864326799153427368592787E-02,
            7.102884842310253397447305465997026228407227220665E-02,
            5.026383572857942403759829860675892897279675661654E-02,
            4.683670010609093810432609684738393586390722052124E-03,
            1.235837891364555000245004813294817451524633100256E-01,
            1.148933497158144016800199601785309838604146040215E-01,
            1.252575774226122633391477702593585307254527198070E-02,
            1.239572396231834242194189674243818619042280816640E-01,
            2.501306413750310579525950767549691151739047969345E-02,
            4.915957918146130094258849161350510503556792927578E-02,
            2.259167374956474713302030584548274729936249753832E-02,
            // weights for Q41
            6.362762978782724559269342300509058175967124446839E-02,
            9.950065827346794643193261975720606296171462239514E-02,
            7.048220002718565366098742295389607994441704889441E-02,
            6.512297339398335645872697307762912795346716454337E-02,
            3.998229150313659724790527138690215186863915308702E-02,
            3.456512257080287509832054272964315588028252136044E-02,
            2.212167975884114432760321569298651047876071264944E-03,
            8.140326425945938045967829319725797511040878579808E-02,
            6.583213447600552906273539578430361199084485578379E-02,
            2.592913726450792546064232192976262988065252032902E-02,
            1.187141856692283347609436153545356484256869129472E-01,
            5.999947605385971985589674757013565610751028128731E-02,
            5.500937980198041736910257988346101839062581489820E-02,
            5.264422421764655969760271538981443718440340270116E-03,
            1.533126874056586959338368742803997744815413565014E-02,
            3.527159369750123100455704702965541866345781113903E-02,
            5.000556431653955124212795201196389006184693561679E-02,
            5.744164831179720106340717579281831675999717767532E-02,
            1.598823797283813438301248206397233634639162043386E-02,
            2.635660410220884993472478832884065450876913559421E-02,
            1.196003937945541091670106760660561117114584656319E-02
    };

    template<class T>
    const size_t Integrator<T>::mw[] = {
            0, 7, 17, 31
    };

    template<class T>
    const size_t Integrator<T>::nw[] = {
            7, 10, 14, 21
    };

}

#endif // INTEGRATOR_H
