//! @file ecom.cxx
//! Program to extract the command line from an Especia result HTML file
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <iostream>

using namespace std;

/// Extracts the command from Especia result HTML. Reads from standard
/// input and writes to standard output.
///
/// @return an exit code.
///
/// @remark Usage: ecom < {result file} [> {target file}]
int main() {
    bool found = false;
    string s;

    while (getline(cin, s)) {
        if (found and s != "</command>") {
            cout << s << endl;
        } else {
            found = (s == "<command>");
        }
    }
    return 0;
}
