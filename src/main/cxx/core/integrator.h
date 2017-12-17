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
    template<class T = real>
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
        explicit Integrator(Formula p = Q27, Formula q = Q41) : p(p), q(q) {
        }

        /**
         * The destructor.
         */
        ~Integrator() = default;

        /**
         * Computes the integral of a function, with the limits supplied as argument, i.e.
         * @f[ \int_{a}^{b} f(x) dx @f].
         *
         * @tparam F The integrand type.
         *
         * @param[in] f The integrand.
         * @param[in] a The lower limit of integration.
         * @param[in] b The upper limit of integration.
         * @param[in] accuracy_goal The (absolute) accuracy goal.
         * @param[in] max_iteration The maximum number of iterations.
         * @return the value of the integral.
         */
        template<class F>
        T integrate(const F &f, T a, T b, T accuracy_goal = T(1.0E-6), natural max_iteration = 100) const {
            Partition<F, Part<F>> partition(f, a, b, p, q);

            for (natural i = 0; i < max_iteration; ++i) {
                if (partition.absolute_error() < accuracy_goal) {
                    break;
                }
                partition.refine();
            }

            return partition.result();
        }

        /**
         * Computes the value of the semi-infinite integral of a function, i.e.
         * @f[ \int_{0}^{\infty} f(x) dx @f].
         *
         * Makes the variable transformation @f$ u = \exp(-x) @f$ and computes
         * @f[ \int_{0}^{1} \frac{f(-\log(u))}{u} du @f].
         *
         * @tparam F The integrand type.
         *
         * @param[in] f The integrand.
         * @param[in] accuracy_goal The (absolute) accuracy goal.
         * @param[in] max_iteration The maximum number of iterations.
         * @return the value of the semi-infinite integral.
         *
         * @attention The integrand must converge rapidly (faster than @f$ 1/x @f$) to zero at infinity
         * @f[ \lim_{x\to\infty} \frac{f(x)}{x} = 0 @f].
         */
        template<class F>
        T integrate_semi_infinite(const F &f, T accuracy_goal = T(1.0E-6), natural max_iteration = 100) const {
            using std::log;

            return integrate([&f](T u) -> T { return u > T(0.0) ? f(-log(u)) / u : T(0.0); }, // infinity maps to zero
                             T(0.0), T(1.0), accuracy_goal, max_iteration);
        }

    private:
        /**
         * A part of a numerical integral.
         *
         * @tparam F The integrand type.
         */
        template<class F>
        class Part {
        public:
            /**
             * The constructor.
             *
             * @param f The integrand.
             * @param a The lower limit of integration.
             * @param b The upper limit of integration.
             * @param p The formula with less quadrature points.
             * @param q The formula with more quadrature points.
             */
            Part(const F &f, T a, T b, Formula p, Formula q)
                    : f(f), a(a), b(b), p(p), q(q), c(T(0.5) * (a + b)), h(T(0.5) * (b - a)), yl(21), yu(21) {
                evaluate();
            }

            /**
             * The destructor.
             */
            ~Part() = default;

            /**
             * Returns the absolute error of the integration result of this part.
             *
             * @return the absolute error of the integration result.
             */
            T absolute_error() const {
                return err;
            }

            /**
             * Returns the integration result of this part.
             *
             * @return the integration result.
             */
            T result() const {
                return res;
            }

            /**
             * Creates a new part from the lower half of this part.
             *
             * @return the lower half part.
             */
            Part *new_lower_part() const {
                auto *part = new Part(this, a, c);

                part->yu[0] = yl[2];
                part->yu[1] = yl[7];
                part->yu[2] = yl[1];
                part->yu[4] = f(part->c + Integrator::xi[4] * part->h);
                part->yu[5] = f(part->c + Integrator::xi[5] * part->h);
                part->yu[6] = yl[0];
                part->yl[0] = yl[2];
                part->yl[1] = yl[8];
                part->yl[2] = yl[3];
                part->yl[3] = yl[4];
                part->yl[4] = yl[5];
                part->yl[5] = yl[9];
                part->yl[6] = yl[6];
                if (nl > 10) {
                    part->yu[3] = yl[10];
                    part->yl[7] = yl[11];
                    part->yl[8] = yl[12];
                    part->yl[9] = yl[13];
                    if (nl > 14) {
                        part->yu[7] = yl[15];
                        part->yu[8] = yl[14];
                        part->yu[9] = f(part->c + Integrator::xi[9] * part->h);
                        part->yu[10] = yl[16];
                        part->yl[10] = yl[17];
                        part->yl[11] = yl[18];
                        part->yl[12] = yl[19];
                        part->yl[13] = yl[20];
                        part->nu = 11;
                        part->nl = 14;
                    } else {
                        part->nu = 7;
                        part->nl = 10;
                    }
                } else {
                    part->yu[3] = f(part->c + Integrator::xi[3] * part->h);
                    part->nu = 7;
                    part->nl = 7;
                }

                part->evaluate();

                return part;
            }

            /**
             * Creates a new part from the upper half of this part.
             *
             * @return the upper half part.
             */
            Part *new_upper_part() const {
                auto *part = new Part(this, c, b);

                part->yl[0] = yu[2];
                part->yl[1] = yu[7];
                part->yl[2] = yu[1];
                part->yl[4] = f(part->c - Integrator::xi[4] * part->h);
                part->yl[5] = f(part->c - Integrator::xi[5] * part->h);
                part->yl[6] = yu[0];
                part->yu[0] = yu[2];
                part->yu[1] = yu[8];
                part->yu[2] = yu[3];
                part->yu[3] = yu[4];
                part->yu[4] = yu[5];
                part->yu[5] = yu[9];
                part->yu[6] = yu[6];
                if (nu > 10) {
                    part->yl[3] = yu[10];
                    part->yu[7] = yu[11];
                    part->yu[8] = yu[12];
                    part->yu[9] = yu[13];
                    if (nu > 14) {
                        part->yl[7] = yu[15];
                        part->yl[8] = yu[14];
                        part->yl[9] = f(part->c - Integrator::xi[9] * part->h);
                        part->yl[10] = yu[16];
                        part->yu[10] = yu[17];
                        part->yu[11] = yu[18];
                        part->yu[12] = yu[19];
                        part->yu[13] = yu[20];
                        part->nl = 11;
                        part->nu = 14;
                    } else {
                        part->nl = 7;
                        part->nu = 10;
                    }
                } else {
                    part->yl[3] = f(part->c - Integrator::xi[3] * part->h);
                    part->nl = 7;
                    part->nu = 7;
                }

                part->evaluate();

                return part;
            }

        private:
            /**
             * Constructs a new part from a parent part.
             *
             * @param parent The parent part.
             * @param a The lower limit of integration.
             * @param b The upper limit of integration.
             */
            Part(const Part *parent, T a, T b)
                    : f(parent->f), a(a), b(b), p(parent->p), q(parent->q),
                      c(T(0.5) * (a + b)), h(T(0.5) * (b - a)), yl(21), yu(21) {
                // do not evaluate
            }

            /**
             * Evaluates the integration result of this part and its absolute error.
             */
            void evaluate() {
                using std::abs;

                res = evaluate(q);
                err = abs(res - evaluate(p));
            }

            /**
             * Evaluates the integration result of this part using the quadrature formula supplied as argument.
             *
             * @param q The quadrature formula.
             * @return the result.
             */
            T evaluate(Formula q) {
                const natural m = Integrator::mw[q];
                const natural n = Integrator::nw[q];

                T result = T(0.0);
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
                for (natural i = 0; i < n; ++i) {
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

            /**
             * The integrand.
             */
            const F &f;

            /**
             * The lower limit of integration.
             */
            const T a;

            /**
             * The upper limit of integration.
             */
            const T b;

            /**
             * The selected formula with less quadrature points.
             */
            const Formula p;

            /**
             * The selected formula with more quadrature points.
             */
            const Formula q;

            /**
             * The center of the interval of integration.
             */
            const T c;

            /**
             * The width of the interval of integration.
             */
            const T h;

            /**
             * The integrand values for the lower half interval of integration.
             */
            std::valarray<T> yl;

            /**
             * The integrand values for the upper half interval of integration.
             */
            std::valarray<T> yu;

            /**
             * The number of evaluated integrand values for the lower half interval.
             */
            natural nl = 0;

            /**
             * The number of evaluated integrand values for the upper half interval.
             */
            natural nu = 0;

            /**
             * The absolute error of the integration result.
             */
            T err = T(0.0);

            /**
             * The integration result.
             */
            T res = T(0.0);
        };

        /**
         * Compares the absolute error of two parts of a numerical integration.
         *
         * @tparam P The part type.
         */
        template<class P>
        class Part_Compare {
        public:
            /**
             * Compares the absolute error of two parts of a numerical integration.
             *
             * @param p The first part.
             * @param q The other part.
             * @return @c true, if the absolute error of the first part is less than that of the other part.
             */
            bool operator()(const P *p, const P *q) const {
                return p->absolute_error() < q->absolute_error();
            }
        };

        /**
         * A partition of a numerical integral into a complete set of disjoint parts.
         *
         * @tparam F The integrand type.
         * @tparam P The part type.
         */
        template<class F, class P>
        class Partition {
        public:
            /**
             * The constructor.
             *
             * @param f The integrand.
             * @param a The lower limit of integration.
             * @param b The upper limit of integration.
             * @param p The formula with less quadrature points.
             * @param q The formula with more quadrature points.
             */
            Partition(const F &f, T a, T b, Formula p, Formula q) : part_compare(Part_Compare<P>()), parts() {
                using std::make_heap;
                using std::push_heap;

                auto *part = new P(f, a, b, p, q);
                make_heap(parts.begin(), parts.end(), part_compare);
                add_part(part);
            }

            /**
             * The destructor.
             */
            ~Partition() { // NOLINT
                for (auto part : parts) {
                    delete part;
                }
            }

            /**
             * Returns the absolute error of the integration result for this partition.
             *
             * @return the absolute error of the integration result.
             */
            T absolute_error() const {
                T err = T(0.0);

                for (auto part : parts) {
                    err += part->absolute_error();
                }

                return err;
            }

            /**
             * Returns the integration result for this partition.
             *
             * @return the integration result.
             */
            T result() const {
                T res = T(0.0);

                for (auto part : parts) {
                    res += part->result();
                }

                return res;
            }

            /**
             * Refines this partition.
             */
            void refine() {
                P *popped = pop_part();
                add_part(popped->new_lower_part());
                add_part(popped->new_upper_part());

                delete popped;
            }

        private:
            /**
             * Removes the part with the largest absolute error of integration from the partition.
             *
             * @return the part with the largest absolute error.
             */
            P *pop_part() {
                using std::pop_heap;

                P *popped = parts.front();
                pop_heap(parts.begin(), parts.end(), part_compare);
                parts.pop_back();

                return popped;
            }

            /**
             * Adds a new part to the partition.
             *
             * @param part The part.
             */
            void add_part(P *part) {
                using std::push_heap;

                parts.push_back(part);
                push_heap(parts.begin(), parts.end(), part_compare);
            }

            /**
             * Compares the absolute error of integration of two parts.
             */
            const Part_Compare<P> part_compare;

            /**
             * The parts of this partition.
             */
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
        static const natural mw[];

        /**
         * The number of quadrature weights.
         */
        static const natural nw[];
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
            T(1.303262173284849021810473057638590518409112513421E-01L),
            T(2.390632866847646220320329836544615917290026806242E-01L),
            T(2.630626354774670227333506083741355715758124943143E-01L),
            T(2.186819313830574175167853094864355208948886875898E-01L),
            T(2.757897646642836865859601197607471574336674206700E-02L),
            T(1.055750100538458443365034879086669791305550493830E-01L),
            T(1.571194260595182254168429283636656908546309467968E-02L),
            // weights for Q19
            T(1.298751627936015783241173611320651866834051160074E-01L),
            T(2.249996826462523640447834514709508786970828213187E-01L),
            T(1.680415725925575286319046726692683040162290325505E-01L),
            T(1.415567675701225879892811622832845252125600939627E-01L),
            T(1.006482260551160175038684459742336605269707889822E-01L),
            T(2.510604860724282479058338820428989444699235030871E-02L),
            T(9.402964360009747110031098328922608224934320397592E-03L),
            T(5.542699233295875168406783695143646338274805359780E-02L),
            T(9.986735247403367525720377847755415293097913496236E-02L),
            T(4.507523056810492466415880450799432587809828791196E-02L),
            // weights for Q27
            T(6.300942249647773931746170540321811473310938661469E-02L),
            T(1.261383225537664703012999637242003647020326905948E-01L),
            T(1.273864433581028272878709981850307363453523117880E-01L),
            T(8.576500414311820514214087864326799153427368592787E-02L),
            T(7.102884842310253397447305465997026228407227220665E-02L),
            T(5.026383572857942403759829860675892897279675661654E-02L),
            T(4.683670010609093810432609684738393586390722052124E-03L),
            T(1.235837891364555000245004813294817451524633100256E-01L),
            T(1.148933497158144016800199601785309838604146040215E-01L),
            T(1.252575774226122633391477702593585307254527198070E-02L),
            T(1.239572396231834242194189674243818619042280816640E-01L),
            T(2.501306413750310579525950767549691151739047969345E-02L),
            T(4.915957918146130094258849161350510503556792927578E-02L),
            T(2.259167374956474713302030584548274729936249753832E-02L),
            // weights for Q41
            T(6.362762978782724559269342300509058175967124446839E-02L),
            T(9.950065827346794643193261975720606296171462239514E-02L),
            T(7.048220002718565366098742295389607994441704889441E-02L),
            T(6.512297339398335645872697307762912795346716454337E-02L),
            T(3.998229150313659724790527138690215186863915308702E-02L),
            T(3.456512257080287509832054272964315588028252136044E-02L),
            T(2.212167975884114432760321569298651047876071264944E-03L),
            T(8.140326425945938045967829319725797511040878579808E-02L),
            T(6.583213447600552906273539578430361199084485578379E-02L),
            T(2.592913726450792546064232192976262988065252032902E-02L),
            T(1.187141856692283347609436153545356484256869129472E-01L),
            T(5.999947605385971985589674757013565610751028128731E-02L),
            T(5.500937980198041736910257988346101839062581489820E-02L),
            T(5.264422421764655969760271538981443718440340270116E-03L),
            T(1.533126874056586959338368742803997744815413565014E-02L),
            T(3.527159369750123100455704702965541866345781113903E-02L),
            T(5.000556431653955124212795201196389006184693561679E-02L),
            T(5.744164831179720106340717579281831675999717767532E-02L),
            T(1.598823797283813438301248206397233634639162043386E-02L),
            T(2.635660410220884993472478832884065450876913559421E-02L),
            T(1.196003937945541091670106760660561117114584656319E-02L)
    };

    template<class T>
    const natural Integrator<T>::mw[] = {
            0, 7, 17, 31
    };

    template<class T>
    const natural Integrator<T>::nw[] = {
            7, 10, 14, 21
    };
}

#endif // INTEGRATOR_H

