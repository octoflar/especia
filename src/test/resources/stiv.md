# Example model definition file
 
Copyright (c) 2016 Ralf Quast

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

    especiv 27182 40 80 0.5 0.000001 2000 10 < stiv.md
    
Based on observations made with the NASA/ESA Hubble Space Telescope, obtained
from the data archive at the Space Telescope Institute. STScI is operated by
the association of Universities for Research in Astronomy, Inc. under the NASA
contract NAS 5-26555.

## Section: sub-damped Lyman-alpha line toward HE 0515-4414

Voigt profiles are used here. The line parameters are:

1. laboratory wavelength (Å)
2. oscillator strength
3. cosmological redshift
4. radial velocity (km s-1)
5. line broadening velocity (km s-1)
6. decadic logarithm of the particle column number density (cm-2)
7. damping constant (s-1)

````````
{
% section 1
% id            source                  begin       end         polynomials
  H_I_1216      stis1216.dat            2596.0      2638.0      3
%
% resolution
% initial       min         max         optimize    reference
  10000         9           11          0
%
% lines
% id
% initial       min         max         optimize    reference

  H_I_1216_01
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           -1800.0     -1400.0     1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_02
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           -1200.0     -800.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_03
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           -800.0     -200.0       1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_DLA
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           -100.0      100.0       1
  0             10.0        200.0       1
  0             19.0        21.0        1
  6.265e+08     0           0           0

  H_I_1216_04
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           200.0       600.0       1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_05
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           1000.0      1400.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_06
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           1300.0      1600.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_07
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           1600.0      1800.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_08
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           1800.0      2000.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_09
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           2000.0      2100.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0

  H_I_1216_10
  1215.6700     0           0           0
  0.4164        0           0           0
  1.1508        0           0           0
  0.0           2200.0      2400.0      1
  0             10.0        200.0       1
  0             12.0        18.0        1
  6.265e+08     0           0           0
}
````````

## References

Reimers, Dieter; Baade, Robert; Quast, Ralf; Levshakov, Sergei A. (2002): *Detection of molecular hydrogen at z = 1.15 toward HE 0515-4414.* 
Astronomy and Astrophysics 410 (3) 785. doi: [10.1051/0004-6361:20031313](http://dx.doi.org/10.1051/0004-6361:20031313).
