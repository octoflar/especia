//! @file emes.cxx
//! Program to extract messages from an Especia result HTML file
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <iostream>

using namespace std;

/// Extracts status messages from Especia result HTML. Reads from standard
/// input and writes to standard output.
///
/// @return an exit code.
///
/// @remark Usage: emes < {result file} [> {target file}]
int main() {
    bool found = false;
    string s;

    while (getline(cin, s)) {
        if (found and not (s == "</message>"
                           // for compatibility with former versions
                           or s == "</mesg>")) {
            cout << s << endl;
        } else {
            found = (s == "<message>"
                     // for compatibility with former versions
                     or s == "<mesg>");
        }
    }

    return 0;
}
