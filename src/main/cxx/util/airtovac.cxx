/// @file airtovac.cxx
/// Utility to convert photon wavelength from air to vacuum
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include <exception>
#include <stdexcept>

#include "../core/dataio.h"
#include "../core/equations.h"
#include "../core/exitcodes.h"

using namespace std;

using especia::natural;
using especia::real;

/// Utility to convert photon wavelength (Angstrom) in spectroscopic data from
/// air to vacuum.
///
/// Further reading:
///
/// B. Edlen (1966).
///  *The refractive index of air.*
///  Metrologia, 2, 2, 71-80.
///  http://dx.doi.org/10.1088/0026-1394/2/2/002
///
/// B. Edlen (1953).
///  *The dispersion of standard air.*
///  Journal of the Optical Society of America, 43, 5, 339.
///
/// @param argc The number of command line arguments supplied.
/// @param argv[0] The program name.
/// @param argv[1] The number of lines to skip at the beginning (optional,
/// default = 0).
/// @return an exit code.
///
/// @remark Usage: airtovac [lines to skip] < {source file} [> {target file}]
int
main (int argc, char *argv[])
{
  using especia::Equations;

  const string program_name (argv[0]);

  try
    {
      if (argc != 1 and argc != 2)
        {
          throw invalid_argument (
              "Error: an invalid number of arguments was supplied");
        }

      natural skip = 0;

      if (argc == 2)
        {
          skip = especia::convert<natural> (string (argv[1]));
        }

      valarray<real> x;
      valarray<real> y;
      valarray<real> z;

      if (especia::get (cin, x, y, z, skip))
        {
          for (size_t i = 0; i < x.size (); ++i)
            {
              x[i] = real (10.0)
                     / especia::solve (Equations::edlen66, real (10.0) / x[i],
                                       real (10.0) / x[i], real (1.0E-08));
            }
          especia::put (cout, x, y, z);
        }
      else
        {
          throw runtime_error ("Error: an input error occurred");
        }
      return 0;
    }
  catch (logic_error &e)
    {
      cerr << e.what () << endl;
      return especia::Exit_Codes::logic_error;
    }
  catch (runtime_error &e)
    {
      cerr << e.what () << endl;
      return especia::Exit_Codes::runtime_error;
    }
  catch (exception &e)
    {
      cerr << e.what () << endl;
      return especia::Exit_Codes::unspecific_exception;
    }
}
