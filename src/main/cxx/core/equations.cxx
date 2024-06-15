/// @file equations.cxx
/// Equations from scientific literature.
/// @author Ralf Quast
/// @date 2024
/// @copyright MIT License
#include "equations.h"

void
especia::Equations::birch94 (const real &x, real &y, real &z)
{
  const real n = 1.0 + 8.34254E-05 + 2.406147E-08 / (130.0E-06 - x * x)
                 + 1.5998E-10 / (38.9E-06 - x * x);
  const real m = (4.812294E-08 * x) / sq (130.0E-06 - x * x)
                 + (3.1996E-10 * x) / sq (38.9E-06 - x * x);

  y = x * n;
  z = n + x * m;
}

void
especia::Equations::edlen53 (const real &x, real &y, real &z)
{
  const real n = 1.0 + 6.43280E-05 + 2.5540E-10 / (0.0000410 - x * x)
                 + 2.949810E-08 / (0.000146 - x * x);
  const real m = (5.1080E-10 * x) / sq (0.0000410 - x * x)
                 + (5.89962E-08 * x) / sq (0.000146 - x * x);

  y = x * n;
  z = n + x * m;
}

void
especia::Equations::edlen66 (const real &x, real &y, real &z)
{
  const real n = 1.0 + 8.34213E-05 + 1.5997E-10 / (0.0000389 - x * x)
                 + 2.406030E-08 / (0.000130 - x * x);
  const real m = (3.1994E-10 * x) / sq (0.0000389 - x * x)
                 + (4.81206E-08 * x) / sq (0.000130 - x * x);

  y = x * n;
  z = n + x * m;
}
