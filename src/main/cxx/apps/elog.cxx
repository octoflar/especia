/// @file elog.cxx
/// Program to extract the log data from an Especia result HTML file
/// @author Ralf Quast
/// @date 2021
/// @copyright MIT License
#include <iostream>

using namespace std;

/// Extracts the log data from Especia result HTML. Reads from standard
/// input and writes to standard output.
///
/// @return an exit code.
///
/// @remark Usage: elog < {result file} [> {target file}]
int main() {
    bool found = false;
    string s;

    while (getline(cin, s)) {
        if (found and s != "</log>") {
            cout << s << endl;
        } else {
            found = (s == "<log>");
        }
    }

    return 0;
}
