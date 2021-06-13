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
    Real evalPoly(const Index degree, const Real x) { return recurrence(degree, x, 1.0, x); }
    
    Real evalPolyDeriv(const Index degree, const Real x) { return recurrence(degree, x, 0.0, 1.0); }

  private:  
    Real recurrence(const Index degree, const Real x, const Real seed0, const Real seed1)
    {
      Real p_n2 = seed0;
      Real p_n1 = seed1;

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
};


#endif
