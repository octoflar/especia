/// @file model.h
/// Parametric model for fitting absorption line regions.
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
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

namespace especia
{

/// The parametric model for fitting absorption line regions.
///
/// @tparam Function The type of profile function.
template <class Function> class Model
{
public:
  /// A bounded constraint.
  ///
  /// @tparam T The number type.
  template <class T = real> class Bounded_Constraint
  {
  public:
    /// Constructs a new strict-bound prior constraint.
    ///
    /// @param[in] lower_bounds The lower bounds.
    /// @param[in] upper_bounds The upper bounds.
    /// @param[in] n The number of bounds.
    Bounded_Constraint (const T lower_bounds[], const T upper_bounds[],
                        natural n)
        : a (lower_bounds, n), b (upper_bounds, n)
    {
    }

    /// The destructor.
    ~Bounded_Constraint () = default;

    /// Tests if a given parameter vector violates the constraint.
    ///
    /// @param[in] x The parameter vector.
    /// @param[in] n The number of parameters to test.
    /// @return @c true, if the parameter vector violates the constraint.
    bool
    is_violated (const T x[], natural n) const
    {
      for (natural i = 0; i < n; ++i)
        {
          if (x[i] < a[i] || x[i] > b[i])
            {
              return true;
            }
        }
      return false;
    }

    /// Computes the cost associated with the constraint.
    ///
    /// @param[in] x The parameter vector.
    /// @param[in] n The number of parameters to take account of.
    /// @return always zero.
    T
    cost (const T x[], natural n) const
    {
      return T (0);
    }

  private:
    const std::valarray<T> a;
    const std::valarray<T> b;
  };

  /// Initializes the model by reading a model definition from in input stream.
  ///
  /// @param is The input stream.
  /// @param os The output stream used for reporting.
  /// @param[in] comment_mark The character to mark a comment.
  /// @param[in] begin_of_section The character to mark the biginning of a
  /// model section.
  /// @param[in] end_of_section The character to mark the end of a model
  /// section.
  /// @return the input stream.
  std::istream &
  get (std::istream &is, std::ostream &os, char comment_mark = '%',
       char begin_of_section = '{', char end_of_section = '}')
  {
    using namespace std;

    typedef map<string, natural>::const_iterator id_index_map_ci;

    const char errmsg[] = "especia::Model<>::get(): Error: ";
    const char dlimsg[] = "duplicate line identifier";
    const char dsimsg[] = "duplicate section identifier";
    const char fnfmsg[] = "file not found";
    const char infmsg[] = "input failed";
    const char srfmsg[] = "self reference";
    const char synmsg[] = "syntax error";
    const char rnfmsg[] = "reference not found";

    vector<especia::Section> sections;

    vector<natural> isc;
    vector<natural> nle;
    vector<natural> nli;

    vector<real> val;
    vector<real> lo;
    vector<real> up;

    vector<bool> msk;
    vector<natural> ind;

    map<string, natural> profile_name_map;
    map<string, natural> section_name_map;

    vector<string> ref;

    stringstream ss;
    stringstream st;
    string line;

    os << "<!DOCTYPE html>\n";
    os << "<html>\n";
    os << "<!--\n";
    os << "<model>\n";

    // Read all lines and write them to the output stream.
    while (readline (is, line))
      {
        ss << line << endl;
        os << line << endl;
      }

    os << "</model>\n";
    os << "-->\n";
    os << "</html>\n";

    // Now strip empty lines and comments
    while (readline (ss, line, comment_mark))
      {
        st << line << '\n';
      }

    natural i = 0;
    size_t j = 0;

    while (getline (st, line, end_of_section))
      if (!st.eof ())
        {
          j = line.find (begin_of_section);

          if (j != string::npos)
            {
              istringstream ist (line.substr (j + 1));

              string sid, pid, fn, s2;
              real a, b;
              natural p;

              // Parse section head
              if (ist >> sid >> fn >> a >> b >> p and getline (ist, s2))
                {
                  if (section_name_map.find (sid) == section_name_map.end ())
                    {
                      section_name_map[sid] = sections.size ();

                      ifstream ifs (fn.c_str ());

                      if (ifs)
                        {
                          especia::Section s;

                          if (s.get (ifs, a, b))
                            {
                              istringstream is2 (s2);
                              while (is2 >> a >> b)
                                s.mask (a, b);

                              sections.push_back (s);
                              isc.push_back (i);
                              nle.push_back (p);
                            }
                          else
                            {
                              is.setstate (ios_base::badbit
                                           | ios_base::failbit);
                              cerr << errmsg << fn << ": " << infmsg << endl;

                              return is;
                            }
                        }
                      else
                        {
                          is.setstate (ios_base::badbit | ios_base::failbit);
                          cerr << errmsg << fn << ": " << fnfmsg << endl;

                          return is;
                        }
                    }
                  else
                    {
                      is.setstate (ios_base::badbit | ios_base::failbit);
                      cerr << errmsg << sid << ": " << dsimsg << endl;

                      return is;
                    }
                }
              else
                {
                  is.setstate (ios_base::badbit | ios_base::failbit);
                  cerr << errmsg << infmsg << endl;

                  return is;
                }

              // Read resolution parameter specification
              if (read (ist, val, lo, up, msk, ref, 1, '\n', true))
                ++i;
              else
                {
                  is.setstate (ios_base::badbit | ios_base::failbit);
                  cerr << errmsg << infmsg << endl;

                  return is;
                }

              natural k = 0;

              // Read profile function parameter specification
              while (ist >> pid)
                if (profile_name_map.find (pid) == profile_name_map.end ())
                  {
                    profile_name_map[pid] = i;

                    if (read (ist, val, lo, up, msk, ref,
                              Function::parameter_count (), '\n', true))
                      {
                        i += Function::parameter_count ();
                        k += 1;
                      }
                    else
                      {
                        is.setstate (ios_base::badbit | ios_base::failbit);
                        cerr << errmsg << infmsg << endl;

                        return is;
                      }
                  }
                else
                  {
                    is.setstate (ios_base::badbit | ios_base::failbit);
                    cerr << errmsg << pid << ": " << dlimsg << endl;

                    return is;
                  }

              nli.push_back (k);
            }
          else
            {
              is.setstate (ios_base::badbit | ios_base::failbit);
              cerr << errmsg << synmsg << endl;

              return is;
            }
        }

    if (!st.bad () and st.eof ())
      {
        // Index independent parameters
        for (natural i = 0, k = 0; i < msk.size (); ++i)
          if (msk[i] and ref[i].empty ())
            {
              if (lo[i] > up[i])
                swap (lo[i], up[i]);

              ind.push_back (k++);
            }
          else
            {
              lo[i] = 0.0;
              up[i] = 0.0;

              ind.push_back (0);
            }

        // Dereference resolution parameter references
        for (id_index_map_ci i = section_name_map.begin ();
             i != section_name_map.end (); ++i)
          {
            const natural j = isc[i->second];

            while (!ref[j].empty ())
              {
                if (section_name_map.find (ref[j]) != section_name_map.end ())
                  {
                    const natural k = isc[section_name_map[ref[j]]];

                    if (j != k)
                      {
                        if (ref[k].empty ())
                          {
                            val[j] = val[k];

                            lo[j] = lo[k];
                            up[j] = up[k];

                            msk[j] = msk[k];
                            ind[j] = ind[k];
                          }
                        ref[j] = ref[k];
                      }
                    else
                      {
                        is.setstate (ios_base::failbit);
                        cerr << errmsg << ref[j] << srfmsg << endl;

                        return is;
                      }
                  }
                else
                  {
                    is.setstate (ios_base::failbit);
                    cerr << errmsg << ref[j] << rnfmsg << endl;

                    return is;
                  }
              }
          }

        // Dereference line parameter references
        for (id_index_map_ci i = profile_name_map.begin ();
             i != profile_name_map.end (); ++i)
          for (natural j = 0; j < Function::parameter_count (); ++j)
            {
              const natural k = i->second + j;

              while (!ref[k].empty ())
                {
                  if (profile_name_map.find (ref[k])
                      != profile_name_map.end ())
                    {
                      const natural l = profile_name_map[ref[k]] + j;

                      if (k != l)
                        {
                          if (ref[l].empty ())
                            {
                              val[k] = val[l];

                              lo[k] = lo[l];
                              up[k] = up[l];

                              msk[k] = msk[l];
                              ind[k] = ind[l];
                            }
                          ref[k] = ref[l];
                        }
                      else
                        {
                          is.setstate (ios_base::failbit);
                          cerr << errmsg << ref[j] << srfmsg << endl;

                          return is;
                        }
                    }
                  else
                    {
                      is.setstate (ios_base::failbit);
                      cerr << errmsg << ref[j] << rnfmsg << endl;

                      return is;
                    }
                }
            }

        const natural m = sections.size ();
        const natural n = msk.size ();

        this->sections = sections;

        this->isc.resize (m);
        this->nle.resize (m);
        this->nli.resize (m);

        copy (isc.begin (), isc.end (), &(this->isc[0]));
        copy (nle.begin (), nle.end (), &(this->nle[0]));
        copy (nli.begin (), nli.end (), &(this->nli[0]));

        this->val.resize (n);
        this->unc.resize (n);

        copy (val.begin (), val.end (), &(this->val[0]));

        this->lo.resize (n);
        this->up.resize (n);

        copy (lo.begin (), lo.end (), &(this->lo[0]));
        copy (up.begin (), up.end (), &(this->up[0]));

        this->msk.resize (n);
        this->ind.resize (n);

        copy (msk.begin (), msk.end (), &(this->msk[0]));
        copy (ind.begin (), ind.end (), &(this->ind[0]));

        this->section_name_map = section_name_map;
        this->profile_name_map = profile_name_map;

        is.clear (is.rdstate () & ~ios_base::failbit);
      }
    else
      is.setstate (ios_base::badbit | ios_base::failbit);

    return is;
  }

  /// Writes the model to an output stream.
  ///
  /// @param os The output stream.
  /// @return the output stream.
  std::ostream &
  put (std::ostream &os) const
  {
    using namespace std;

    typedef map<string, natural>::const_iterator id_index_map_ci;

    const ios_base::fmtflags fmt = os.flags ();

    os.setf (ios_base::fmtflags ());
    os.setf (ios_base::fixed, ios_base::floatfield);
    os.setf (ios_base::left, ios_base::adjustfield);

    os << "<!DOCTYPE html>\n";
    os << "<html>\n";
    os << "<!--\n";
    os << "<data>\n";
    os << sections;
    os << "</data>\n";
    os << "-->\n";

    os << "<head>\n";
    os << "  <title>Parameter Table</title>\n";
    os << "</head>\n";
    os << "<body>\n";
    os << "<table border=\"1\" cellspacing=\"2\" cellpadding=\"2\" "
          "width=\"100%\">\n";
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

    for (auto i = section_name_map.begin (); i != section_name_map.end (); ++i)
      {
        const natural j = i->second;

        const string id = i->first;
        const size_t px = sections[j].valid_data_count ();
        const real st = sections[j].cost ();

        os.precision (2);

        os << "    <tr>\n";
        os << "      <td>" << id << "</td>\n";
        os << "      <td>" << sections[j].lower_bound () << "</td>\n";
        os << "      <td>" << sections[j].upper_bound () << "</td>\n";
        os << "      <td>" << nle[j] << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 2, isc[j]);
        os << "</td>\n";
        os << "      <td>" << px << "</td>\n";
        os << "      <td><strong>" << st << "</strong></td>\n";
        os << "      <td>" << st / px << "</td>\n";
        os << "    </tr>\n";
      }

    os << "  </tbody>\n";
    os << "</table>\n";
    os << "<br>\n";
    os << "<table border=\"1\" cellspacing=\"2\" cellpadding=\"2\" "
          "width=\"100%\">\n";
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
    if (Function::parameter_count () == 8)
      {
        os << "      <td>&Delta;&alpha;/&alpha;<br>(10<sup>-6</sup>)</td>\n";
      }
    os << "    </tr>\n";
    os << "  </thead>\n";
    os << "  <tbody align=\"left\">\n";

    const Equivalent_Width_Calculator<Integrator<real> > calculator;

    for (auto i = profile_name_map.begin (); i != profile_name_map.end (); ++i)
      {
        const natural j = i->second;
        const string id = i->first;

        const real c = 1.0E-3 * speed_of_light;
        const real x = val[j];
        const real z = val[j + 2];
        const real v = val[j + 3];
        const real w = x * (1.0 + z) * (1.0 + v / c);
        const real dx = unc[j];
        const real dz = unc[j + 2];
        const real dv = unc[j + 3];
        const real dw
            = dx
              + x * sqrt (sq ((1.0 + v / c) * dz) + sq ((1.0 + z) * dv / c));
        const real ew = calculator.calculate (Function (&val[j]), milli);

        os.precision (4);

        os << "    <tr>\n";
        os << "      <td>" << id << "</td>\n";
        os << "      <td>" << w << " &plusmn; " << dw << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 4, j);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::scientific, 3, j + 1);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 7, j + 2);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 3, j + 3);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 3, j + 4);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 3, j + 5);
        os << "</td>\n";
        os << "      <td>";
        put_parameter (os, ios_base::fixed, 3, ew);
        os << "</td>\n";
        if (Function::parameter_count () == 8)
          {
            os << "      <td>";
            put_parameter (os, ios_base::fixed, 3, j + 7);
            os << "</td>\n";
          }
        os << "    </tr>\n";
      }

    os << "  </tbody>\n";
    os << "</table>\n";

    os << "<address>" << endl;
    os << " Created by <cite>" << project_title << "</cite>. "
       << project_doi_html << "<br>" << endl;
    os << " " << project_long_name << "<br>" << endl;
    os << " " << system_name << "<br>" << endl;
    os << " " << cxx_compiler << " " << cxx_compiler_version << "<br>" << endl;
    os << "</address>" << endl;

    os << "</body>\n";
    os << "</html>\n";

    os.flush ();
    os.flags (fmt);

    return os;
  }

  /// Computes the value of the cost function for given (extrinsic) model
  /// parameter values.
  ///
  /// @param[in] x The (extrinsic) model parameter values.
  /// @param[in] n The number of (extrinsic) model parameters.
  /// @return the value of the cost function.
  real
  operator() (const real x[], natural n) const
  {
    return cost (x, n);
  }

  /// Sets new model parameters values and uncertainties.
  ///
  /// @param[in] x The (extrinsic) model parameter values.
  /// @param[in] u The (extrinsic) model parameter uncertainties.
  void
  set (const real x[], const real u[])
  {
    for (natural i = 0; i < val.size (); ++i)
      {
        if (msk[i])
          {
            val[i] = x[ind[i]];
            unc[i] = u[ind[i]];
          }
        else
          {
            unc[i] = 0.0;
          }
      }
    for (natural i = 0; i < sections.size (); ++i)
      {
        sections[i].apply (nle[i], val[isc[i]],
                           Superposition<Function> (nli[i], &val[isc[i] + 1]));
      }
  }

  /// Computes the value of the cost function for given (extrinsic) model
  /// parameter values.
  ///
  /// @param[in] x The (extrinsic) model parameter values.
  /// @param[in] n The number of (extrinsic) model parameters.
  /// @return the value of the cost function.
  real
  cost (const real x[], natural n) const
  {
    using std::valarray;

    valarray<real> y = val;
    for (natural i = 0; i < y.size (); ++i)
      {
        if (msk[i])
          {
            y[i] = x[ind[i]];
          }
      }
    real d = 0.0;
    for (natural i = 0; i < sections.size (); ++i)
      {
        d += sections[i].cost (
            Superposition<Function> (nli[i], &y[isc[i] + 1]), y[isc[i]],
            nle[i]);
      }
    return d;
  }

  /// Returns the number of (extrinsic) model parameters.
  ///
  /// @return the number of (extrinsic) model parameters.
  natural
  get_parameter_count () const
  {
    return ind.max () + 1;
  }

  /// Returns the initial (extrinsic) parameter values.
  ///
  /// @return the initial (extrinsic) parameter values.
  std::valarray<real>
  get_initial_parameter_values () const
  {
    std::valarray<real> x (get_parameter_count ());

    for (natural i = 0, j = 0; i < msk.size (); ++i)
      {
        if (msk[i] and ind[i] == j)
          {
            x[j++] = 0.5 * (lo[i] + up[i]);
          }
      }

    return x;
  }

  /// Returns the initial step sizes associated with (extrinsic) parameter
  /// values.
  ///
  /// @return the initial step sizes associated with (extrinsic) parameter
  /// values.
  std::valarray<real>
  get_initial_local_step_sizes () const
  {
    std::valarray<real> z (get_parameter_count ());

    for (natural i = 0, j = 0; i < msk.size (); ++i)
      {
        if (msk[i] and ind[i] == j)
          {
            z[j++] = 0.5 * (up[i] - lo[i]);
          }
      }

    return z;
  }

  /// Returns the bound constraints associated with (extrinsic) model
  /// parameters.
  ///
  /// @return the bound constraints associated with (extrinsic) model
  /// parameters.
  Bounded_Constraint<real>
  get_constraint () const
  {
    std::valarray<real> a (get_parameter_count ());
    std::valarray<real> b (get_parameter_count ());

    for (natural i = 0, j = 0; i < msk.size (); ++i)
      {
        if (msk[i] and ind[i] == j)
          {
            a[j] = lo[i];
            b[j] = up[i];
            ++j;
          }
      }

    return Bounded_Constraint<real> (&a[0], &b[0], get_parameter_count ());
  }

private:
  std::ostream &
  put_parameter (std::ostream &os, std::ios_base::fmtflags f, natural p,
                 real parameter) const
  {
    using namespace std;

    const ios_base::fmtflags fmt = os.flags ();

    os.setf (f, ios_base::floatfield);
    os.precision (p);

    os << parameter;

    os.flags (fmt);

    return os;
  }

  std::ostream &
  put_parameter (std::ostream &os, std::ios_base::fmtflags f, natural p,
                 natural parameter_index) const
  {
    using namespace std;

    const ios_base::fmtflags fmt = os.flags ();

    os.setf (f, ios_base::floatfield);
    os.precision (p);

    os << val[parameter_index];
    if (msk[parameter_index])
      os << " &plusmn; " << unc[parameter_index];

    os.flags (fmt);

    return os;
  }

  /// The container of model sections.
  std::vector<especia::Section> sections;

  /// The index of the first model parameter of each section.
  std::valarray<natural> isc;
  /// The polynomial degree of the continuum approximation in each section.
  std::valarray<natural> nle;
  /// The number of spectral lines in each section.
  std::valarray<natural> nli;

  /// The intrinsic model parameters.
  std::valarray<real> val;
  /// The uncertainty associated with intrinsic model parameters.
  std::valarray<real> unc;
  /// The lower bounds of intrinsic model parameters.
  std::valarray<real> lo;
  /// The upper bounds of intrinsic model parameters.
  std::valarray<real> up;

  /// The mask is @c true when a parameter is optimized, @false when it is
  /// considered constant.
  std::valarray<bool> msk;
  /// The index maps from extrinsic to intrinsic model parameterss.
  std::valarray<natural> ind;

  std::map<std::string, natural> section_name_map;
  std::map<std::string, natural> profile_name_map;
};

}

#endif // ESPECIA_MODEL_H
