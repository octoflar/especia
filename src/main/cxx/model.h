// Parametric model for fitting absorption line regions
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
#ifndef RQ_MODEL_H
#define RQ_MODEL_H

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include "optimize.h"
#include "profiles.h"
#include "readline.h"
#include "section.h"


namespace RQ {
    template<class profile_function>
    class model;
}

template<class profile_function>
class RQ::model {
public:
    std::istream &get(std::istream &is, std::ostream &os, char comment_mark = '%', char begin_of_section = '{',
                      char end_of_section = '}');

    std::ostream &put(std::ostream &os) const;

    void compute_model(const double x[], size_t n);

    double statistics(const double x[], size_t n) const;

    template<class normal_deviate, class sym_eig_decomp>
    bool optimize(size_t parent_number,
                  size_t population_size,
                  double step_size,
                  const double accuracy_goal,
                  unsigned long stop_generation,
                  unsigned trace,
                  normal_deviate &ndev, sym_eig_decomp &eig, std::ostream &os);

private:
    std::ostream &put_parameter(std::ostream &os, std::ios_base::fmtflags f,
                                int precision, size_t parameter_index) const;

    std::vector<RQ::section> sec;

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


template<class profile_function>
std::istream &
RQ::model<profile_function>::get(std::istream &is, std::ostream &os, char cm, char bos, char eos) {
    using namespace std;

    typedef map<string, size_t>::const_iterator id_index_map_ci;

    const char errmsg[] = "RQ::model<>::get(): Error: ";
    const char dlimsg[] = "duplicate line identifier";
    const char dsimsg[] = "duplicate section identifier";
    const char fnfmsg[] = "file not found";
    const char infmsg[] = "input failed";
    const char srfmsg[] = "self reference";
    const char synmsg[] = "syntax error";
    const char rnfmsg[] = "reference not found";

    vector<RQ::section> sec;

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
    string s;

    os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
    os << "<html>\n";
    os << "<!--\n";
    os << "<model>\n";

    // Strip empty lines. Write non-empty lines to standard output.
    while (readline(is, s)) {
        ss << s << '\n';
        os << s << '\n';
    }

    os << "</model>\n";
    os << "-->\n";
    os << "</html>\n";

    // Now strip comments
    while (readline(ss, s, cm)) {
        st << s << '\n';
    }

    size_t i = 0;
    size_t j = 0;

    while (getline(st, s, eos))
        if (!st.eof()) {
            j = s.find(bos);

            if (j != string::npos) {
                istringstream ist(s.substr(j + 1));

                string sid, pid, fn, s2;
                double a, b;
                size_t p;

                // Parse section head
                if (ist >> sid >> fn >> a >> b >> p and getline(ist, s2)) {
                    if (sim.find(sid) == sim.end()) {
                        sim[sid] = sec.size();

                        ifstream ifs(fn.c_str());

                        if (ifs) {
                            RQ::section s;

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

                size_t j = 0;

                // Read profile function parameter specification
                while (ist >> pid)
                    if (pim.find(pid) == pim.end()) {
                        pim[pid] = i;

                        if (read(ist, val, lo, up, msk, ref, profile_function::parameters, '\n', true)) {
                            i += profile_function::parameters;
                            j += 1;
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

                nli.push_back(j);
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
            for (size_t j = 0; j < profile_function::parameters; ++j) {
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

template<class profile_function>
std::ostream &
RQ::model<profile_function>::put(std::ostream &os) const {
    using namespace std;

    typedef map<string, size_t>::const_iterator id_index_map_ci;

    const ios_base::fmtflags fmt = os.flags();

    os.setf(ios_base::fmtflags());
    os.setf(ios_base::fixed, ios_base::floatfield);
    os.setf(ios_base::left, ios_base::adjustfield);

    os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
    os << "<html>\n";
    os << "<!--\n";
    os << "  Created by\n";
    os << "  Evolutionary spectrum inversion and analysis (Especia)\n";
    os << "  http://dx.doi.org/10.6084/m9.figshare.4167999\n";
    os << "-->\n";
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
    os << "      <td>Legendre<br>Polynomials</td>\n";
    os << "      <td>Resolution</td>\n";
    os << "      <td>Data Points</td>\n";
    os << "      <td>Cost</td>\n";
    os << "      <td>Cost per<br>Data Point</td>\n";
    os << "    </tr>\n";
    os << "  </thead>\n";
    os << "  <tbody align=\"left\">\n";

    for (id_index_map_ci i = sim.begin(); i != sim.end(); ++i) {
        const size_t j = i->second;

        const string id = i->first;
        const size_t px = sec[j].selection_size();
        const double st = sec[j].cost();

        os.precision(2);

        os << "    <tr>\n";
        os << "      <td>" << id << "</td>\n";
        os << "      <td>" << sec[j].begin() << "</td>\n";
        os << "      <td>" << sec[j].end() << "</td>\n";
        os << "      <td>" << nle[j] << "</td>\n";
        os << "      <td>";
        put_parameter(os, ios_base::fixed, 2, isc[j]);
        os << "</td>\n";
        os << "      <td>" << px << "</td>\n";
        os << "      <td>" << st << "</td>\n";
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
#if defined(RQ_MANY_MULTIPLET_ANALYSIS)
    os << "      <td>&Delta;&alpha/&alpha;<br>(10<sup>-5</sup>)</td>\n";
#endif
    os << "    </tr>\n";
    os << "  </thead>\n";
    os << "  <tbody align=\"left\">\n";

    for (id_index_map_ci i = pim.begin(); i != pim.end(); ++i) {
        const size_t j = i->second;

        const string id = i->first;

        const double x = val[j];
        const double z = val[j + 2];
        const double v = val[j + 3];
        const double w = x * (1.0 + z) * (1.0 + v / 299792.458);
        const double dz = err[j + 2];
        const double dv = err[j + 3];
        const double dw = x * sqrt(RQ::sqr((1.0 + v / 299792.458) * dz) + RQ::sqr((1.0 + z) * dv / 299792.458));

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
#if defined(RQ_MANY_MULTIPLET_ANALYSIS)
        os << "      <td>";
        put_parameter(os, ios_base::fixed, 3, j + 7);
        os << "</td>\n";
#endif
        os << "    </tr>\n";
    }

    os << "  </tbody>\n";
    os << "</table>\n";
    os << "</body>\n";
    os << "</html>\n";

    os.flush();
    os.flags(fmt);

    return os;
}

template<class profile_function>
std::ostream &
RQ::model<profile_function>::put_parameter(std::ostream &os, std::ios_base::fmtflags f,
                                           int p, size_t i) const {
    using namespace std;

    const ios_base::fmtflags fmt = os.flags();

    os.setf(f, ios_base::floatfield);
    os.precision(p);

    os << val[i];
    if (msk[i])
        os << " &plusmn; " << err[i];

    os.flags(fmt);

    return os;
}

template<class profile_function>
void
RQ::model<profile_function>::compute_model(const double x[], size_t n) {
    for (size_t i = 0; i < val.size(); ++i)
        if (msk[i])
            val[i] = x[ind[i]];

    for (size_t i = 0; i < sec.size(); ++i)
        sec[i].compute_model(RQ::superposition<profile_function>(nli[i], &val[isc[i] + 1]),
                             nle[i], val[isc[i]]);
}

template<class profile_function>
double
RQ::model<profile_function>::statistics(const double x[], size_t n) const {
    using namespace std;

    valarray<double> y = val;
    for (size_t i = 0; i < y.size(); ++i)
        if (msk[i])
            y[i] = x[ind[i]];

    double d = 0.0;

    for (size_t i = 0; i < sec.size(); ++i)
        d += sec[i].cost(RQ::superposition<profile_function>(nli[i], &y[isc[i] + 1]),
                         nle[i], y[isc[i]]);

    return d;
}

template<class profile_function>
template<class normal_deviate, class sym_eig_decomp>
bool
RQ::model<profile_function>::optimize(size_t parent_number,
                                      size_t population_size,
                                      double step_size,
                                      const double accuracy_goal,
                                      unsigned long stop_generation,
                                      unsigned trace,
                                      normal_deviate &ndev, sym_eig_decomp &eig, std::ostream &os) {
    using namespace std;

    const char beglog[] = "<log>";
    const char endlog[] = "</log>";
    const char begmsg[] = "<message>";
    const char endmsg[] = "</message>";
    const char optmsg[] = "RQ::model<>::optimize(): Message: Optimization completed";
    const char stpmsg[] = "RQ::model<>::optimize(): Warning: Optimization stopped at generation ";
    const char uflmsg[] = "RQ::model<>::optimize(): Warning: Mutation variance underflow";

    valarray<double> w(1.0, parent_number);
    for (size_t i = 0; i < parent_number; ++i)
        w[i] = log((parent_number + 1.0) / (i + 1));

    const size_t n = ind.max() + 1;
    const double var_parent_number = RQ::sqr(w.sum()) / w.apply(RQ::sqr).sum();
    const double cs = (var_parent_number + 2.0) / (var_parent_number + n + 3.0);
    const double cc = 4.0 / (n + 4.0);
    const double acov = 1.0 / var_parent_number;
    const double ccov = acov * (2.0 / RQ::sqr(n + sqrt(2.0))) + (1.0 - acov) * min(1.0,
                                                                                   (2.0 * var_parent_number - 1.0) /
                                                                                   (RQ::sqr(n + 2.0) +
                                                                                    var_parent_number));
    const double step_size_damping = cs + 1.0 + 2.0 * max(0.0, sqrt((var_parent_number - 1.0) / (n + 1.0)) - 1.0);

    valarray<double> pc(0.0, n);
    valarray<double> ps(0.0, n);

    valarray<double> d(0.0, n);
    valarray<double> B(0.0, RQ::sqr(n));
    valarray<double> C = B;

    valarray<double> x(n);
    valarray<double> a(n);
    valarray<double> b(n);

    for (size_t i = 0, j = 0; i < msk.size(); ++i)
        if (msk[i] and ind[i] == j) {
            a[j] = lo[i];
            b[j] = up[i];
            x[j] = 0.5 * (a[j] + b[j]);
            ++j;
        }

    for (size_t i = 0, ii = 0; i < n; ++i, ii += n + 1) {
        d[i] = 0.5 * (b[i] - a[i]);
        B[ii] = 1.0;
        C[ii] = d[i] * d[i];
    }

    bool is_opt = false;
    bool is_ufl = false;

    unsigned long g = 0;
    double z = 0.0;

    os << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
    os << "<html>\n";
    os << "<!--\n";
    os << "  Created by\n";
    os << "  Evolutionary spectrum inversion and analysis (Especia)\n";
    os << "  http://dx.doi.org/10.6084/m9.figshare.4167999\n";
    os << "-->\n";

    if (trace > 0) {
        os << "<!--" << endl;
        os << beglog << endl;
    }

    while (!is_opt and !is_ufl and g < stop_generation) {
        unsigned long stop = stop_generation;

        if (trace > 0)
            stop = min(g + trace, stop_generation);

        RQ::optimize(this, &model<profile_function>::statistics, &x[0], n,
                     &a[0], &b[0],
                     parent_number,
                     population_size,
                     &w[0],
                     step_size,
                     step_size_damping,
                     cs,
                     cc,
                     ccov,
                     acov,
                     accuracy_goal,
                     stop,
                     &d[0],
                     &B[0],
                     &C[0],
                     &ps[0],
                     &pc[0],
                     g,
                     is_opt,
                     is_ufl,
                     z,
                     ndev, eig, less<double>());

        if (trace > 0) {
            const ios_base::fmtflags fmt = os.flags();

            const int p = 4;
            const int w = 12;

            os.setf(ios_base::fmtflags());
            os.setf(ios_base::scientific, ios_base::floatfield);
            os.setf(ios_base::right, ios_base::adjustfield);
            os.precision(p);

            os << setw(8) << g;
            os << setw(w) << z;
            os << setw(w) << step_size * d.min();
            os << setw(w) << step_size * d.max();
            os << endl;

            os.flags(fmt);
        }
    }

    if (trace > 0) {
        os << endlog << endl;
        os << "-->" << endl;
    }

    os << "<!--" << endl;
    os << begmsg << endl;

    if (is_opt)
        os << optmsg << endl;
    else if (is_ufl)
        os << uflmsg << endl;
    else
        os << stpmsg << g << endl;

    os << endmsg << endl;
    os << "-->" << endl;
    os << "</html>\n";

    // Compute uncertainty
    if (is_opt)
        RQ::scale_step_size(this, &model<profile_function>::statistics, &x[0], n, step_size, &d[0], &B[0]);
    for (size_t i = 0, ii = 0; i < n; ++i, ii += n + 1)
        d[i] = step_size * sqrt(C[ii]);
    for (size_t i = 0; i < msk.size(); ++i)
        if (msk[i])
            err[i] = d[ind[i]];
        else
            err[i] = 0.0;

    compute_model(&x[0], n);

    return (is_opt or is_ufl or g >= stop_generation);
}

#endif // RQ_MODEL_H
