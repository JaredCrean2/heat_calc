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
    Real recurrence(const Index degree, const Real x);

    Real recurrence_dx(const Index degree, const Real x);
};


#endif
