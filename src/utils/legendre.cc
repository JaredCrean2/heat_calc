#include "utils/legendre.h"

Real LegendrePoly::recurrence(const Index degree, const Real x)
{
  Real p_n2 = 1.0;
  Real p_n1 = x;

  if (degree == 0)
    return p_n2;

  if (degree == 1)
    return p_n1;

  Real p_n;
  for (int n=2; n <= degree; ++n)
  {
    p_n = ((2*n - 1)*x*p_n1 - (n - 1)*p_n2)/n;
    p_n2 = p_n1;
    p_n1 = p_n;
  }

  return p_n;
}

Real LegendrePoly::recurrence_dx(const Index degree, const Real x)
{
  Real p_n2 = 1.0, p_n2_dx = 0.0;
  Real p_n1 = x,   p_n1_dx = 1.0;

  if (degree == 0)
    return p_n2_dx;

  if (degree == 1)
    return p_n1_dx;

  Real p_n, p_n_dx;
  for (int n=2; n <= degree; ++n)
  {
    p_n    = ((2*n - 1)*x*p_n1 - (n - 1)*p_n2)/n;
    p_n_dx = ((2*n - 1)*(p_n1 + x*p_n1_dx) - (n - 1)*p_n2_dx)/n;

    p_n2 = p_n1; p_n2_dx = p_n1_dx;
    p_n1 = p_n;  p_n1_dx = p_n_dx;
  }

  return p_n_dx;
}

