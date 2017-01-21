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

    especid 27182 10 20 0.5 0.000001 2000 10 < uvesd.md
    
Data based on observations made with ESO Telescopes at the La Silla
or Paranal Observatories under programme ID 066.A-0212.

## Section: fine-structure multiplet of CI 1560

Doppler profiles are used here. The line parameters are:

1. laboratory wavelength (Å)
2. oscillator strength
3. cosmological redshift
4. radial velocity (km s-1)
5. line broadening velocity (km s-1)
6. decadic logarithm of the particle column number density (cm-2)

````````
{
% section 1
% id            source                  begin       end         polynomials
  C_I_1560      uves1561.dat            3355.00     3359.00     3
%
% spectral resolution
% initial       min         max         optimize    reference
  50000         40000       60000       1
%
% absorption lines
% id
% initial       min         max         optimize    reference
  C_I_3P0-3D1_1 % C I
  1560.3092     0           0           0
  0.0719        0           0           0
  1.1508        1.150       1.152       0
  0             -10.0       0.0         1
  0             0.0         10.0        1
  12.0          11.0        14.0        1
  
  C_I_3P1-3D2_1 % C I*
  1560.6822     0           0           0
  0.0539        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  12.0          11.0        14.0        1
  
  C_I_3P1-3D1_1 % C I*
  1560.7090     0           0           0
  0.0180        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P1-3D2_1
  
  C_I_3P2-3D2_1 % C I**
  1561.3402     0           0           0
  0.0108        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  12.0          11.0        14.0        1
  
  C_I_3P2-3D1_1 % C I**
  1561.3667     0           0           0
  0.000716      0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P2-3D2_1
  
  C_I_3P2-3D3_1 % C I**
  1561.4384     0           0           0
  0.0603        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P2-3D2_1
  
  C_I_3P0-3D1_2 % C I
  1560.3092     0           0           0
  0.0719        0           0           0
  1.1508        1.150       1.152       0
  0             -10.0       10.0        1
  0             0.0         10.0        1
  12.0          11.0        14.0        1
  
  C_I_3P1-3D2_2 % C I*
  1560.6822     0           0           0
  0.0539        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  12.0          11.0        14.0        1
  
  C_I_3P1-3D1_2 % C I*
  1560.7090     0           0           0
  0.0180        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P1-3D2_2
  
  C_I_3P2-3D2_2 % C I**
  1561.3402     0           0           0
  0.0108        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  12.0          11.0        14.0        1
  
  C_I_3P2-3D1_2 % C I**
  1561.3667     0           0           0
  0.000716      0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P2-3D2_2
  
  C_I_3P2-3D3_2 % C I**
  1561.4384     0           0           0
  0.0603        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P2-3D2_2
}
````````

## Section: fine-structure multiplet of CI 1656

This section exhibits three unidentified lines that affect the
modelling of the background continuum, if not considered. All
parameters of the CI 1656 lines are tied to those of CI 1560.

````````
{
% section 2
% id            source                  begin       end         polynomials
  C_I_1656      uves1657.dat            3562.00     3567.00     3
%
% spectral resolution
% initial       min         max         optimize    reference
  0             0           0           0           C_I_1560
%
% absorption lines
% id
% initial       min         max         optimize    reference
  C_I_3P1-3P2_1 % C I*
  1656.2672     0           0           0
  0.0589        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P1-3D2_1
  
  C_I_3P0-3P1_1 % C I
  1656.9283     0           0           0
  0.139         0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  
  C_I_3P2-3P2_1 % C I**
  1657.0082     0           0           0
  0.104         0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P2-3D2_1
  
  C_I_3P1-3P1_1 % C I*
  1657.3792     0           0           0
  0.0356        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P1-3D2_1
  
  C_I_3P1-3P0_1 % C I*
  1657.9068     0           0           0
  0.0473        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P1-3D2_1
  
  C_I_3P2-3P1_1 % C I**
  1658.1212     0           0           0
  0.0356        0           0           0
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P0-3D1_1
  0             0           0           0           C_I_3P2-3D2_1
  
  C_I_3P1-3P2_2 % C I*
  1656.2672     0           0           0
  0.0589        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P1-3D2_2
  
  C_I_3P0-3P1_2 % C I
  1656.9283     0           0           0
  0.139         0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  
  C_I_3P2-3P2_2 % C I**
  1657.0082     0           0           0
  0.104         0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P2-3D2_2
  
  C_I_3P1-3P1_2 % C I*
  1657.3792     0           0           0
  0.0356        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P1-3D2_2
  
  C_I_3P1-3P0_2 % C I*
  1657.9068     0           0           0
  0.0473        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P1-3D2_2
  
  C_I_3P2-3P1_2 % C I**
  1658.1212     0           0           0
  0.0356        0           0           0
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P0-3D1_2
  0             0           0           0           C_I_3P2-3D2_2
  
  G1 % generic
  3564.10       3563.90     3564.20     1
  1.0           0           0           0
  0.0           0           0           0
  0.0           0           0           0
  7.0           4.0         10.0        1
  12.0          10.0        14.0        1
  
  G2 % generic
  3564.30       3564.20     3564.40     1
  1.0           0           0           0
  0.0           0           0           0
  0.0           0           0           0
  0             0           0           0           G1
  12.0          10.0        14.0        1
  
  G3 % generic
  3564.55       3564.40     3564.70     1
  1.0           0           0           0
  0.0           0           0           0
  0.0           0           0           0
  0             0           0           0           G1
  12.0          10.0        14.0        1
}
````````

## References

Quast, Ralf; Baade, Robert; Reimers, Dieter (2002): *Fine-structure diagnostics of neutral carbon toward HE 0515-4414.* Astronomy and Astrophysics 386 (3) 796. doi: [10.1051/0004-6361:20020342](http://dx.doi.org/10.1051/0004-6361:20020342).

