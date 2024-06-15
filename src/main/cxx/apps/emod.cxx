/// @file emod.cxx
/// Program to extract the model definition from an Especia result HTML file
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include <iostream>

using namespace std;

/// Extracts the model definition from Especia result HTML. Reads from standard
/// input and writes to standard output.
///
/// @return an exit code.
///
/// @remark Usage: emod < {result file} [> {target file}]
int
main ()
{
  bool found = false;
  string s;

  while (getline (cin, s))
    {
      if (found
          and not(s == "</model>"
                  // for compatibility with former versions
                  or s == "</job>" or s == "</input>"))
        {
          cout << s << endl;
        }
      else
        {
          found = (s == "<model>"
                   // for compatibility with former versions
                   or s == "<job>" or s == "<input>");
        }
    }

  return 0;
}
