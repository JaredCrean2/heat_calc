#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "ProjectDefs.h"

// class for evaluating the Legendre polynomials
// Uses Bonnets recursion formula: n P_n(x) = (2n - 1) P_{n-1}(x) - (n - 1) P_{n-2}(x)
// where n is a zero-based index
class LegendrePoly
{
  using Index = int;

  public:
    Real evalPoly(const Index degree, const Real x) { return recurrence(degree, x); }
    
    Real evalPolyDeriv(const Index degree, const Real x) { return recurrence_dx(degree, x); }

  private:  
    Real recurrence(const Index degree, const Real x)
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

    Real recurrence_dx(const Index degree, const Real x)
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
};


#endif
