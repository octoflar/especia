# Example model definition file

Copyright (c) 2018 Ralf Quast

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## About

To optimize this model type, for instance:

    especid 27182 4 8 0.5 0.0001 1000 10 < na.md

## Section: Na-D doublet

Doppler profiles are used here. The line parameters are:

1. laboratory wavelength (Å)
2. oscillator strength
3. cosmological redshift
4. radial velocity (km s-1)
5. line broadening velocity (km s-1)
6. decadic logarithm of the particle column number density (cm-2)

Consult the Especia tutorial on [model definition files](https://github.com/octoflar/especia/wiki/Model-definition-files)
for further explanation.

```
{
% section 1
% id            source                  begin       end         polynomials
  Na_D          resources/spec000.dat   5860.0      5910.0      2
%
% spectral resolution (1E+3)
% initial       min         max         optimize    reference
  3             2           6           1
%
% absorption lines
% id
% initial       min         max         optimize    reference    comment
  Na_D_1
  5889.95       5889.0      5891.0      1
  0.01          0           0           0
  0             0           0           0
  0             0           0           0
  0             2.0         80.0        1
  0             12.0        15.0        1

  Na_D_2
  5895.92       5895.0      5897.0      1
  0.01          0           0           0
  0             0           0           0
  0             0           0           0
  0             2.0         80.0        1
  0             12.0        15.0        1
}
```

## Thanks

Based on a test case provided by Martin Wendt, AIP, Germany.

