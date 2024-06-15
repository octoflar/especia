/// @file readline.cxx
/// Data input procedures.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include <cctype>

#include "readline.h"

using namespace std;

istream &
especia::readline (istream &is, string &s, char comment_mark, char eol)
{
  bool empty_line = true;

  while (empty_line and getline (is, s, eol) and comment_mark != '\0')
    {
      size_t i = 0;

      while (i < s.length ()
             and (comment_mark == '\0' or s[i] != comment_mark))
        empty_line = isspace (s[i++]) and empty_line;

      s.erase (i);
    }

  return is;
}
