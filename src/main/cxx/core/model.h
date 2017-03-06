/// @file model.h
/// Parametric model for fitting absorption line regions.
/// Copyright (c) 2016 Ralf Quast
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
#include "profiles.h"
#include "readline.h"
#include "section.h"
#include "../optimization/optimizer.h"


namespace especia {

    template<class Profile>
    class Model {
    public:

        std::istream &get(std::istream &is,
                          std::ostream &os,
                          char comment_mark = '%',
                          char begin_of_section = '{',
                          char end_of_section = '}') {
            using namespace std;

            typedef map<string, size_t>::const_iterator id_index_map_ci;

            const char errmsg[] = "especia::Model<>::get(): Error: ";
            const char dlimsg[] = "duplicate line identifier";
            const char dsimsg[] = "duplicate section identifier";
            const char fnfmsg[] = "file not found";
            const char infmsg[] = "input failed";
            const char srfmsg[] = "self reference";
            const char synmsg[] = "syntax error";
            const char rnfmsg[] = "reference not found";

            vector<especia::Section> sec;

            vector<size_t> isc;
            vector<size_t> nle;
            vector<size_t> nli;

            vector<double> val;
            vector<double> lo;
            vector<double> up;

            vector<bool> msk;
            vector<size_t> ind;

            map<string, size_t> pim;
            map<string, size_t> sim;

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

            size_t i = 0;
            size_t j = 0;

            while (getline(st, line, end_of_section))
                if (!st.eof()) {
                    j = line.find(begin_of_section);

                    if (j != string::npos) {
                        istringstream ist(line.substr(j + 1));

                        string sid, pid, fn, s2;
                        double a, b;
                        size_t p;

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

                        size_t k = 0;

                        // Read profile function parameter specification
                        while (ist >> pid)
                            if (pim.find(pid) == pim.end()) {
                                pim[pid] = i;

                                if (read(ist, val, lo, up, msk, ref, Profile::parameter_count, '\n', true)) {
                                    i += Profile::parameter_count;
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
                for (size_t i = 0, k = 0; i < msk.size(); ++i)
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
                    const size_t j = isc[i->second];

                    while (!ref[j].empty()) {
                        if (sim.find(ref[j]) != sim.end()) {
                            const size_t k = isc[sim[ref[j]]];

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
                    for (size_t j = 0; j < Profile::parameter_count; ++j) {
                        const size_t k = i->second + j;

                        while (!ref[k].empty()) {
                            if (pim.find(ref[k]) != pim.end()) {
                                const size_t l = pim[ref[k]] + j;

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

                const size_t m = sec.size();
                const size_t n = msk.size();

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

            typedef map<string, size_t>::const_iterator id_index_map_ci;

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
                const size_t j = i->second;

                const string id = i->first;
                const size_t px = sec[j].valid_data_count();
                const double st = sec[j].cost();

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
#if defined(ESPECIA_MANY_MULTIPLET_ANALYSIS)
            os << "      <td>&Delta;&alpha;/&alpha;<br>(10<sup>-6</sup>)</td>\n";
#endif
            os << "    </tr>\n";
            os << "  </thead>\n";
            os << "  <tbody align=\"left\">\n";

            for (id_index_map_ci i = pim.begin(); i != pim.end(); ++i) {
                const size_t j = i->second;
                const string id = i->first;

                const double c = 1.0E-3 * speed_of_light;
                const double x = val[j];
                const double z = val[j + 2];
                const double v = val[j + 3];
                const double w = x * (1.0 + z) * (1.0 + v / c);
                const double dx = err[j];
                const double dz = err[j + 2];
                const double dv = err[j + 3];
                const double dw = dx + x * sqrt(sqr((1.0 + v / c) * dz) + sqr((1.0 + z) * dv / c));

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
            os << " " << PROJECT_LONG_NAME << " " << "<a href=\"" << DOI << "\">" << DOI << "</a>" << "<br>\n";
            os << " " << SYSTEM << " " << "<br>\n";
            os << " " << CXX_COMPILER << " " << CXX_COMPILER_VERSION << "<br>\n";
            os << "</address>\n";
            os << "</body>\n";
            os << "</html>\n";

            os.flush();
            os.flags(fmt);

            return os;
        }

        double operator()(const double x[], size_t n) const {
            return cost(x, n);
        }

        void apply(const double x[], const double z[]) {
            for (size_t i = 0; i < val.size(); ++i) {
                if (msk[i]) {
                    val[i] = x[ind[i]];
                    err[i] = z[ind[i]];
                } else {
                    err[i] = 0.0;
                }
            }
            for (size_t i = 0; i < sec.size(); ++i) {
                sec[i].apply(Superposition<Profile>(nli[i], &val[isc[i] + 1]), val[isc[i]], nle[i]);
            }
        }

        double cost(const double x[], size_t n) const {
            using std::valarray;

            valarray<double> y = val;
            for (size_t i = 0; i < y.size(); ++i) {
                if (msk[i]) {
                    y[i] = x[ind[i]];
                }
            }
            double d = 0.0;
            for (size_t i = 0; i < sec.size(); ++i) {
                d += sec[i].cost(Superposition<Profile>(nli[i], &y[isc[i] + 1]), y[isc[i]], nle[i]);
            }
            return d;
        }

        size_t get_parameter_count() const {
            return ind.max() + 1;
        }

        std::valarray<double> get_initial_parameter_values() const {
            std::valarray<double> x(get_parameter_count());

            for (size_t i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    x[j++] = 0.5 * (lo[i] + up[i]);
                }
            }

            return x;
        }

        std::valarray<double> get_initial_local_step_sizes() const {
            std::valarray<double> z(get_parameter_count());

            for (size_t i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    z[j++] = 0.5 * (up[i] - lo[i]);
                }
            }

            return z;
        }

        Bounded_Constraint<double> get_constraint() const {
            std::valarray<double> a(get_parameter_count());
            std::valarray<double> b(get_parameter_count());

            for (size_t i = 0, j = 0; i < msk.size(); ++i) {
                if (msk[i] and ind[i] == j) {
                    a[j] = lo[i];
                    b[j] = up[i];
                    ++j;
                }
            }

            return Bounded_Constraint<double>(&a[0], &b[0], get_parameter_count());
        }

    private:
        std::ostream &put_parameter(std::ostream &os, std::ios_base::fmtflags f, int p, size_t parameter_index) const {
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

        std::valarray<size_t> isc;
        std::valarray<size_t> nle;
        std::valarray<size_t> nli;

        std::valarray<double> val;
        std::valarray<double> err;
        std::valarray<double> lo;
        std::valarray<double> up;

        std::valarray<bool> msk;
        std::valarray<size_t> ind;

        std::map<std::string, size_t> sim;
        std::map<std::string, size_t> pim;
    };

}

#endif // ESPECIA_MODEL_H
