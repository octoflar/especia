/// @file model.h
/// Parametric model for fitting absorption line regions.
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
#ifndef ESPECIA_MODEL_H
#define ESPECIA_MODEL_H

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

#include "config.h"
#include "integrator.h"
#include "profiles.h"
#include "readline.h"
#include "section.h"

namespace especia {

    template<class Profile>
    class Model {
    public:

        /**
         * A bounded constraint.
         *
         * @tparam T The number type.
         */
        template<class T = Real>
        class Bounded_Constraint {
        public:
            /**
             * Constructs a new strict-bound prior constraint.
             *
             * @param[in] lower_bounds The lower bounds.
             * @param[in] upper_bounds The upper bounds.
             * @param[in] n The number of bounds.
             */
            Bounded_Constraint(const T lower_bounds[], const T upper_bounds[], Natural n)
                    : a(lower_bounds, n), b(upper_bounds, n) {
            }

            /**
             * The destructor.
             */
            ~Bounded_Constraint() {
            }

            /**
             * Tests if a given parameter vector violates the constraint.
             *
             * @param[in] x The parameter vector.
             * @param[in] n The number of parameters to test.
             * @return @c true, if the parameter vector violates the constraint.
             */
            bool is_violated(const T x[], Natural n) const {
                for (Natural i = 0; i < n; ++i) {
                    if (x[i] < a[i] || x[i] > b[i]) {
                        return true;
                    }
                }
                return false;
            }

            /**
             * Computes the cost associated with the constraint.
             *
             * @param[in] x The parameter vector.
             * @param[in] n The number of parameters to take account of.
             * @return always zero.
             */
            T cost(const T x[], Natural n) const {
                return T(0);
            }

        private:
            const std::valarray<T> a;
            const std::valarray<T> b;
        };

        std::istream &get(std::istream &is,
                          std::ostream &os,
                          char comment_mark = '%',
                          char begin_of_section = '{',
                          char end_of_section = '}') {
            using namespace std;

            typedef map<string, Natural>::const_iterator id_index_map_ci;

            const char errmsg[] = "especia::Model<>::get(): Error: ";
            const char dlimsg[] = "duplicate line identifier";
            const char dsimsg[] = "duplicate section identifier";
            const char fnfmsg[] = "file not found";
            const char infmsg[] = "input failed";
            const char srfmsg[] = "self reference";
            const char synmsg[] = "syntax error";
            const char rnfmsg[] = "reference not found";

            vector<especia::Section> sec;

            vector<Natural> isc;
            vector<Natural> nle;
            vector<Natural> nli;

            vector<Real> val;
            vector<Real> lo;
            vector<Real> up;

            vector<bool> msk;
            vector<Natural> ind;

            map<string, Natural> pim;
            map<string, Natural> sim;

            vector<string> ref;

            stringstream ss;
            stringstream st;
            string line;

            os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
            os << "<html>\n";
            os << "<!--\n";
            os << "<model>\n";

            // Read all lines and write them to standard output.
            while (readline(is, line)) {
                ss << line << endl;
                os << line << endl;
            }

            os << "</model>\n";
            os << "-->\n";
            os << "</html>\n";

            // Now strip empty lines and comments
            while (readline(ss, line, comment_mark)) {
                st << line << '\n';
            }

            Natural i = 0;
            size_t j = 0;

            while (getline(st, line, end_of_section))
                if (!st.eof()) {
                    j = line.find(begin_of_section);

                    if (j != string::npos) {
                        istringstream ist(line.substr(j + 1));

                        string sid, pid, fn, s2;
                        Real a, b;
                        Natural p;

                        // Parse section head
                        if (ist >> sid >> fn >> a >> b >> p and getline(ist, s2)) {
                            if (sim.find(sid) == sim.end()) {
                                sim[sid] = sec.size();

                                ifstream ifs(fn.c_str());

                                if (ifs) {
                                    especia::Section s;

                                    if (s.get(ifs, a, b)) {
                                        istringstream is2(s2);
                                        while (is2 >> a >> b)
                                            s.mask(a, b);

                                        sec.push_back(s);
                                        isc.push_back(i);
                                        nle.push_back(p);
                                    } else {
                                        is.setstate(ios_base::badbit | ios_base::failbit);
                                        cerr << errmsg << fn << ": " << infmsg << endl;

                                        return is;
                                    }
                                } else {
                                    is.setstate(ios_base::badbit | ios_base::failbit);
                                    cerr << errmsg << fn << ": " << fnfmsg << endl;

                                    return is;
                                }
                            } else {
                                is.setstate(ios_base::badbit | ios_base::failbit);
                                cerr << errmsg << sid << ": " << dsimsg << endl;

                                return is;
                            }
                        } else {
                            is.setstate(ios_base::badbit | ios_base::failbit);
                            cerr << errmsg << infmsg << endl;

                            return is;
                        }

                        // Read resolution parameter specification
                        if (read(ist, val, lo, up, msk, ref, 1, '\n', true))
                            ++i;
                        else {
                            is.setstate(ios_base::badbit | ios_base::failbit);
                            cerr << errmsg << infmsg << endl;

                            return is;
                        }

                        Natural k = 0;

                        // Read profile function parameter specification
                        while (ist >> pid)
                            if (pim.find(pid) == pim.end()) {
                                pim[pid] = i;

                                if (read(ist, val, lo, up, msk, ref, Profile::get_parameter_count(), '\n', true)) {
                                    i += Profile::get_parameter_count();
                                    k += 1;
                                } else {
                                    is.setstate(ios_base::badbit | ios_base::failbit);
                                    cerr << errmsg << infmsg << endl;

                                    return is;
                                }
                            } else {
                                is.setstate(ios_base::badbit | ios_base::failbit);
                                cerr << errmsg << pid << ": " << dlimsg << endl;

                                return is;
                            }

                        nli.push_back(k);
                    } else {
                        is.setstate(ios_base::badbit | ios_base::failbit);
                        cerr << errmsg << synmsg << endl;

                        return is;
                    }
                }

            if (!st.bad() and st.eof()) {
                // Index independent parameters
                for (Natural i = 0, k = 0; i < msk.size(); ++i)
                    if (msk[i] and ref[i].empty()) {
                        if (lo[i] > up[i])
                            swap(lo[i], up[i]);

                        ind.push_back(k++);
                    } else {
                        lo[i] = 0.0;
                        up[i] = 0.0;

                        ind.push_back(0);
                    }

                // Dereference resolution parameter references
                for (id_index_map_ci i = sim.begin(); i != sim.end(); ++i) {
                    const Natural j = isc[i->second];

                    while (!ref[j].empty()) {
                        if (sim.find(ref[j]) != sim.end()) {
                            const Natural k = isc[sim[ref[j]]];

                            if (j != k) {
                                if (ref[k].empty()) {
                                    val[j] = val[k];

                                    lo[j] = lo[k];
                                    up[j] = up[k];

                                    msk[j] = msk[k];
                                    ind[j] = ind[k];
                                }
                                ref[j] = ref[k];
                            } else {
                                is.setstate(ios_base::failbit);
                                cerr << errmsg << ref[j] << srfmsg << endl;

                                return is;
                            }
                        } else {
                            is.setstate(ios_base::failbit);
                            cerr << errmsg << ref[j] << rnfmsg << endl;

                            return is;
                        }
                    }
                }

                // Dereference line parameter references
                for (id_index_map_ci i = pim.begin(); i != pim.end(); ++i)
                    for (Natural j = 0; j < Profile::get_parameter_count(); ++j) {
                        const Natural k = i->second + j;

                        while (!ref[k].empty()) {
                            if (pim.find(ref[k]) != pim.end()) {
                                const Natural l = pim[ref[k]] + j;

                                if (k != l) {
                                    if (ref[l].empty()) {
                                        val[k] = val[l];

                                        lo[k] = lo[l];
                                        up[k] = up[l];

                                        msk[k] = msk[l];
                                        ind[k] = ind[l];
                                    }
                                    ref[k] = ref[l];
                                } else {
                                    is.setstate(ios_base::failbit);
                                    cerr << errmsg << ref[j] << srfmsg << endl;

                                    return is;
                                }
                            } else {
                                is.setstate(ios_base::failbit);
                                cerr << errmsg << ref[j] << rnfmsg << endl;

                                return is;
                            }
                        }
                    }

                const Natural m = sec.size();
                const Natural n = msk.size();

                this->sec = sec;

                this->isc.resize(m);
                this->nle.resize(m);
                this->nli.resize(m);

                copy(isc.begin(), isc.end(), &(this->isc[0]));
                copy(nle.begin(), nle.end(), &(this->nle[0]));
                copy(nli.begin(), nli.end(), &(this->nli[0]));

                this->val.resize(n);
                this->err.resize(n);

                copy(val.begin(), val.end(), &(this->val[0]));

                this->lo.resize(n);
                this->up.resize(n);

                copy(lo.begin(), lo.end(), &(this->lo[0]));
                copy(up.begin(), up.end(), &(this->up[0]));

                this->msk.resize(n);
                this->ind.resize(n);

                copy(msk.begin(), msk.end(), &(this->msk[0]));
                copy(ind.begin(), ind.end(), &(this->ind[0]));

                this->sim = sim;
                this->pim = pim;

                is.clear(is.rdstate() & ~ios_base::failbit);
            } else
                is.setstate(ios_base::badbit | ios_base::failbit);

            return is;
        }

        std::ostream &put(std::ostream &os) const {
            using namespace std;

            typedef map<string, Natural>::const_iterator id_index_map_ci;

            const ios_base::fmtflags fmt = os.flags();

            os.setf(ios_base::fmtflags());
            os.setf(ios_base::fixed, ios_base::floatfield);
            os.setf(ios_base::left, ios_base::adjustfield);

            os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
            os << "<html>\n";
            os << "<!--\n";
            os << "<data>\n";
            os << sec;
            os << "</data>\n";
            os << "-->\n";

            os << "<head>\n";
            os << "  <title>Parameter Table</title>\n";
            os << "</head>\n";
            os << "<body>\n";
            os << "<table border=\"1\" cellspacing=\"2\" cellpadding=\"2\" width=\"100%\">\n";
            os << "  <thead align=\"center\" valign=\"middle\">\n";
            os << "    <tr>\n";
            os << "      <td>Section</td>\n";
            os << "      <td>Start<br>Wavelength<br>(&Aring;)</td>\n";
            os << "      <td>End<br>Wavelength<br>(&Aring;)</td>\n";
            os << "      <td>Legendre Basis<br>Polynomials</td>\n";
            os << "      <td>Resolution<br>(10<sup>3</sup>)</td>\n";
            os << "      <td>Data Points</td>\n";
            os << "      <td>Cost</td>\n";
            os << "      <td>Cost per<br>Data Point</td>\n";
            os << "    </tr>\n";
            os << "  </thead>\n";
            os << "  <tbody align=\"left\">\n";

            for (id_index_map_ci i = sim.begin(); i != sim.end(); ++i) {
                const Natural j = i->second;

                const string id = i->first;
                const size_t px = sec[j].valid_data_count();
                const Real st = sec[j].cost();

                os.precision(2);

                os << "    <tr>\n";
                os << "      <td>" << id << "</td>\n";
                os << "      <td>" << sec[j].lower_bound() << "</td>\n";
                os << "      <td>" << sec[j].upper_bound() << "</td>\n";
                os << "      <td>" << nle[j] << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 2, isc[j]);
                os << "</td>\n";
                os << "      <td>" << px << "</td>\n";
                os << "      <td><strong>" << st << "</strong></td>\n";
                os << "      <td>" << st / px << "</td>\n";
                os << "    </tr>\n";
            }

            os << "  </tbody>\n";
            os << "</table>\n";
            os << "<br>\n";
            os << "<table border=\"1\" cellspacing=\"2\" cellpadding=\"2\" width=\"100%\">\n";
            os << "  <thead align=\"center\" valign=\"middle\">\n";
            os << "    <tr>\n";
            os << "      <td>Line</td>\n";
            os << "      <td>Observed<br>Wavelength<br>(&Aring;)</td>\n";
            os << "      <td>Rest<br>Wavelength<br>(&Aring;)</td>\n";
            os << "      <td>Oscillator<br>Strength</td>\n";
            os << "      <td>Redshift</td>\n";
            os << "      <td>Radial<br>Velocity<br>(km s<sup>-1</sup>)</td>\n";
            os << "      <td>Broadening<br>Velocity<br>(km s<sup>-1</sup>)</td>\n";
            os << "      <td>Log. Column<br>Density<br>(cm<sup>-2</sup>)</td>\n";
            os << "      <td>Equivalent<br>Width<br>(m&Aring;)</td>\n";
#if defined(ESPECIA_MANY_MULTIPLET_ANALYSIS)
            os << "      <td>&Delta;&alpha;/&alpha;<br>(10<sup>-6</sup>)</td>\n";
#endif
            os << "    </tr>\n";
            os << "  </thead>\n";
            os << "  <tbody align=\"left\">\n";

            const Equivalent_Width_Calculator<Integrator<Real>> calculator;

            for (id_index_map_ci i = pim.begin(); i != pim.end(); ++i) {
                const Natural j = i->second;
                const string id = i->first;

                const Real c = 1.0E-3 * speed_of_light;
                const Real x = val[j];
                const Real z = val[j + 2];
                const Real v = val[j + 3];
                const Real w = x * (1.0 + z) * (1.0 + v / c);
                const Real dx = err[j];
                const Real dz = err[j + 2];
                const Real dv = err[j + 3];
                const Real dw = dx + x * sqrt(sq((1.0 + v / c) * dz) + sq((1.0 + z) * dv / c));
                const Real ew = calculator.calculate(Profile(&val[j]));

                os.precision(4);

                os << "    <tr>\n";
                os << "      <td>" << id << "</td>\n";
                os << "      <td>" << w << " &plusmn; " << dw << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 4, j);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::scientific, 3, j + 1);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 7, j + 2);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 3, j + 3);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 3, j + 4);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 3, j + 5);
                os << "</td>\n";
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 1, ew);
                os << "</td>\n";
#if defined(ESPECIA_MANY_MULTIPLET_ANALYSIS)
                os << "      <td>";
                put_parameter(os, ios_base::fixed, 3, j + 7);
                os << "</td>\n";
#endif
                os << "    </tr>\n";
            }

            os << "  </tbody>\n";
            os << "</table>\n";
            os << "<address>\n";
            os << " Created by Evolutionary spectrum inversion and analysis (Especia).<br>\n";
            os << " " << project_long_name << " " << "<a href=\"" << project_doi << "\">" << project_doi << "</a>" << "<br>\n";
            os << " " << system_name << " " << "<br>\n";
            os << " " << cxx_compiler << " " << cxx_compiler_version << "<br>\n";
            os << "</address>\n";
            os << "</body>\n";
            os << "</html>\n";

            os.flush();
            os.flags(fmt);

            return os;
        }

        Real operator()(const Real x[], Natural n) const {
            return cost(x, n);
        }

        void set(const Real x[], const Real z[]) {
            for (Natural i = 0; i < val.size(); ++i) {
                if (msk[i]) {
                    val[i] = x[ind[i]];
                    err[i] = z[ind[i]];
                } else {
                    err[i] = 0.0;
                }
            }
            for (Natural i = 0; i < sec.size(); ++i) {
                sec[i].apply(Superposition<Profile>(nli[i], &val[isc[i] + 1]), val[isc[i]], nle[i]);
            }
        }

        Real cost(const Real x[], Natural n) const {
            using std::valarray;

            valarray<Real> y = val;
            for (Natural i = 0; i < y.size(); ++i) {
                if (msk[i]) {
                    y[i] = x[ind[i]];
                }
            }
            Real d = 0.0;
            for (Natural i = 0; i < sec.size(); ++i) {
                d += sec[i].cost(Superposition<Profile>(nli[i], &y[isc[i] + 1]), y[isc[i]], nle[i]);
            }
            return d;
        }

        Natural get_parameter_count() const {
            return ind.max() + 1;
        }

        std::valarray<Real> get_initial_parameter_values() const {
            std::valarray<Real> x(get_parameter_count());

            for (Natural i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    x[j++] = 0.5 * (lo[i] + up[i]);
                }
            }

            return x;
        }

        std::valarray<Real> get_initial_local_step_sizes() const {
            std::valarray<Real> z(get_parameter_count());

            for (Natural i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    z[j++] = 0.5 * (up[i] - lo[i]);
                }
            }

            return z;
        }

        Bounded_Constraint<Real> get_constraint() const {
            std::valarray<Real> a(get_parameter_count());
            std::valarray<Real> b(get_parameter_count());

            for (Natural i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    a[j] = lo[i];
                    b[j] = up[i];
                    ++j;
                }
            }

            return Bounded_Constraint<Real>(&a[0], &b[0], get_parameter_count());
        }

    private:
        std::ostream &put_parameter(std::ostream &os, std::ios_base::fmtflags f, Natural p, Real parameter) const {
            using namespace std;

            const ios_base::fmtflags fmt = os.flags();

            os.setf(f, ios_base::floatfield);
            os.precision(p);

            os << parameter;

            os.flags(fmt);

            return os;
        }

        std::ostream &put_parameter(std::ostream &os, std::ios_base::fmtflags f, Natural p, Natural parameter_index) const {
            using namespace std;

            const ios_base::fmtflags fmt = os.flags();

            os.setf(f, ios_base::floatfield);
            os.precision(p);

            os << val[parameter_index];
            if (msk[parameter_index])
                os << " &plusmn; " << err[parameter_index];

            os.flags(fmt);

            return os;
        }

        std::vector<especia::Section> sec;

        std::valarray<Natural> isc;
        std::valarray<Natural> nle;
        std::valarray<Natural> nli;

        std::valarray<Real> val;
        std::valarray<Real> err;
        std::valarray<Real> lo;
        std::valarray<Real> up;

        std::valarray<bool> msk;
        std::valarray<Natural> ind;

        std::map<std::string, Natural> sim;
        std::map<std::string, Natural> pim;
    };

}

#endif // ESPECIA_MODEL_H
