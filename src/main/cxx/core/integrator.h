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

#include <algorithm>
#include <cmath>
#include <limits>
#include <valarray>
#include <vector>

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
         * Constructs a new integrator, based on the formulas supplied as argument.
         *
         * @param[in] p The formula with less quadrature points.
         * @param[in] q The formula with more quadrature points.
         */
        Integrator(Formula p = Q13, Formula q = Q19) : p(p), q(q) {
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
         * @param[in] accuracy_goal The (absolute) accuracy goal.
         * @param[in] max_iteration The maximum number of iterations.
         * @return the value of the integral.
         */
        template<class F>
        T integrate(const F &f, T a, T b, T accuracy_goal = T(1.0E-6), size_t max_iteration = 100) const {
            Partition<F, Part<F>> partition(f, a, b, p, q);

            for (size_t i = 0; i < max_iteration; ++i) {
                if (partition.get_absolute_error() < accuracy_goal) {
                    break;
                }
                partition.refine();
            }

            return partition.get_result();
        }

        /**
         * Computes the value of the semi-infinite integral of a function, i.e.
         * @f[ \int_{0}^{\infty} f(x) dx @f].
         *
         * Makes the variable transformation @f[ x = exp(-u) @f] and computes
         * @f[ \int_{0}^{1} \frac{f(-log(u))}{u} du @f].
         *
         * @tparam F The function type.
         *
         * @param[in] f The function. To be integrable @c f must vanish at infinity.
         * @param[in] accuracy_goal The (absolute) accuracy goal.
         * @param[in] max_iteration The maximum number of iterations.
         * @return the value of the semi-infinite integral.
         */
        template<class F>
        T integrate_semi_infinite(const F &f, T accuracy_goal = T(1.0E-6), size_t max_iteration = 100) const {
            using std::log;
            using std::numeric_limits;

            return integrate([&f](T u) -> T {
                return f(-log(u)) / u; }, numeric_limits<T>::epsilon(), T(1.0), accuracy_goal, max_iteration);
        }

    private:
        /**
         * A part of a numerical integral.
         *
         * @tparam F The function type.
         */
        template<class F>
        class Part {
        public:
            Part(const F &f, T a, T b, Formula p, Formula q)
                    : f(f), a(a), b(b), p(p), q(q), c(T(0.5) * (a + b)), h(T(0.5) * (b - a)), yl(21), yu(21) {
                evaluate();
            }

            ~Part() {

            }

            T get_absolute_error() const {
                return absolute_error;
            }

            T get_result() const {
                return result;
            }

            Part *lower_half() const {
                Part *half = new Part(this, a, c);

                half->yu[0] = yl[2];
                half->yu[1] = yl[7];
                half->yu[2] = yl[1];
                half->yu[4] = f(half->c + Integrator::xi[4] * half->h);
                half->yu[5] = f(half->c + Integrator::xi[5] * half->h);
                half->yu[6] = yl[0];
                half->yl[0] = yl[2];
                half->yl[1] = yl[8];
                half->yl[2] = yl[3];
                half->yl[3] = yl[4];
                half->yl[4] = yl[5];
                half->yl[5] = yl[9];
                half->yl[6] = yl[6];
                if (nl > 10) {
                    half->yu[3] = yl[10];
                    half->yl[7] = yl[11];
                    half->yl[8] = yl[12];
                    half->yl[9] = yl[13];
                    if (nl > 14) {
                        half->yu[7] = yl[15];
                        half->yu[8] = yl[14];
                        half->yu[9] = f(half->c + Integrator::xi[9] * half->h);
                        half->yu[10] = yl[16];
                        half->yl[10] = yl[17];
                        half->yl[11] = yl[18];
                        half->yl[12] = yl[19];
                        half->yl[13] = yl[20];
                        half->nu = 11;
                        half->nl = 14;
                    } else {
                        half->nu = 7;
                        half->nl = 10;
                    }
                } else {
                    half->yu[3] = f(half->c + Integrator::xi[3] * half->h);
                    half->nu = 7;
                    half->nl = 7;
                }

                half->evaluate();

                return half;
            }

            Part *upper_half() const {
                Part *half = new Part(this, c, b);

                half->yl[0] = yu[2];
                half->yl[1] = yu[7];
                half->yl[2] = yu[1];
                half->yl[4] = f(half->c - Integrator::xi[4] * half->h);
                half->yl[5] = f(half->c - Integrator::xi[5] * half->h);
                half->yl[6] = yu[0];
                half->yu[0] = yu[2];
                half->yu[1] = yu[8];
                half->yu[2] = yu[3];
                half->yu[3] = yu[4];
                half->yu[4] = yu[5];
                half->yu[5] = yu[9];
                half->yu[6] = yu[6];
                if (nu > 10) {
                    half->yl[3] = yu[10];
                    half->yu[7] = yu[11];
                    half->yu[8] = yu[12];
                    half->yu[9] = yu[13];
                    if (nu > 14) {
                        half->yl[7] = yu[15];
                        half->yl[8] = yu[14];
                        half->yl[9] = f(half->c - Integrator::xi[9] * half->h);
                        half->yl[10] = yu[16];
                        half->yu[10] = yu[17];
                        half->yu[11] = yu[18];
                        half->yu[12] = yu[19];
                        half->yu[13] = yu[20];
                        half->nl = 11;
                        half->nu = 14;
                    } else {
                        half->nl = 7;
                        half->nu = 10;
                    }
                } else {
                    half->yl[3] = f(half->c - Integrator::xi[3] * half->h);
                    half->nl = 7;
                    half->nu = 7;
                }

                half->evaluate();

                return half;
            }

        private:
            Part(const Part *parent, T a, T b)
                    : f(parent->f), a(a), b(b), p(parent->p), q(parent->q),
                      c(T(0.5) * (a + b)), h(T(0.5) * (b - a)), yl(21), yu(21) {
                // do not evaluate
            }

            void evaluate() {
                using std::abs;

                result = evaluate(q);
                absolute_error = abs(result - evaluate(p));
            }

            T evaluate(Formula q) {
                const size_t m = Integrator::mw[q];
                const size_t n = Integrator::nw[q];

                T result = T(0.0);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
                for (size_t i = 0; i < n; ++i) {
                    if (i >= nl) {
                        yl[i] = f(c - h * Integrator::xi[i]);
                    }
                    if (i >= nu) {
                        yu[i] = f(c + h * Integrator::xi[i]);
                    }
                    result += (yl[i] + yu[i]) * Integrator::wi[m + i];
                }
                if (nl < n) {
                    nl = n;
                }
                if (nu < n) {
                    nu = n;
                }

                return result * h;
            }

            const F &f;
            const T a;
            const T b;
            const Formula p;
            const Formula q;

            const T c;
            const T h;

            std::valarray<T> yl;
            std::valarray<T> yu;

            size_t nl = 0;
            size_t nu = 0;

            T absolute_error = T(0.0);
            T result = T(0.0);
        };

        /**
         * Compares two integral parts.
         *
         * @tparam P The part type.
         */
        template<class P>
        class Part_Compare {
        public:
            /**
             * Compares two integral parts.
             *
             * @param p The first part.
             * @param q The other part.
             * @return @c true, if the absolute error of the first part is less than that of the other part.
             */
            bool operator()(const P *p, const P *q) const {
                return p->get_absolute_error() < q->get_absolute_error();
            }
        };

        /**
         * The partion of a numerical integral.
         *
         * @tparam F The function type.
         * @tparam P The part type.
         */
        template<class F, class P>
        class Partition {
        public:
            Partition(const F &f, T a, T b, Formula p, Formula q) : part_compare(Part_Compare<P>()), parts() {
                using std::make_heap;
                using std::push_heap;

                P *part = new P(f, a, b, p, q);
                make_heap(parts.begin(), parts.end(), part_compare);
                push_part(part);
            }

            ~Partition() {
                for (auto part : parts) {
                    delete part;
                }
            }

            T get_absolute_error() const {
                T absolute_error = T(0.0);

                for (auto part : parts) {
                    absolute_error += part->get_absolute_error();
                }

                return absolute_error;
            }

            T get_result() const {
                T result = T(0.0);

                for (auto part : parts) {
                    result += part->get_result();
                }

                return result;
            }

            void refine() {
                P *popped = pop_part();
                push_part(popped->lower_half());
                push_part(popped->upper_half());

                delete popped;
            }

        private:
            P *pop_part() {
                using std::pop_heap;

                P *popped = parts.front();
                pop_heap(parts.begin(), parts.end(), part_compare);
                parts.pop_back();

                return popped;
            }

            void push_part(P *part) {
                using std::push_heap;

                parts.push_back(part);
                push_heap(parts.begin(), parts.end(), part_compare);
            }

            const Part_Compare<P> part_compare;

            std::vector<P *> parts;
        };

        /**
         * The selected quadrature formula with less points.
         */
        const Formula p;

        /**
         * The selected quadrature formula with more points.
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
