#include "utils/lagrange.h"

LagrangeBasis::LagrangeBasis(const std::vector<double>& pts) :
  m_pts(pts), m_denoms(pts.size(), 1)
{
  for (Index i=0; i < pts.size(); ++i)
    for (Index j=0; j < pts.size(); ++j)
    {
      if (i == j)
        continue;

      m_denoms[i] *= pts[i] - pts[j];
    }
}


double LagrangeBasis::evalPoly(Index p, double x)
{
  double val = 1;
  for (Index i=0; i < getNumPoints(); ++i)
  {
    if (i == p)
      continue;

    val *= (x - m_pts[i]);
  }

  return val/m_denoms[p];
}


double LagrangeBasis::evalPolyDeriv(Index p, double x)
{
  double val = 0;
  for (Index k=0; k < getNumPoints();  ++k)
  {
    if (k == p)
      continue;

    double val_k = 1;
    for (Index i=0; i < getNumPoints(); ++i)
    {
      if (i == p || i == k)
        continue;

      val_k *= (x - m_pts[i]);
    }

    val += val_k;
  }

  return val/m_denoms[p];
}


