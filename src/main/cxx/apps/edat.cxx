//! @file edat.cxx
//! Program to extract the spectroscopic data from an Especia result HTML file
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <iostream>

using namespace std;

/// Extracts the section data from Especia result HTML. Reads from standard
/// input and writes to standard output.
///
/// @return an exit code.
///
/// @remark Usage: edat < {result file} [> {target file}]
int main() {
    bool found = false;
    string s;

    while (getline(cin, s)) {
        if (found and s != "</data>") {
            cout << s << endl;
        } else {
            found = (s == "<data>");
        }
    }

    return 0;
}
