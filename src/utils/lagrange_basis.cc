#include "utils/lagrange_basis.h"

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



int LagrangeMemoizer::getIdx(const std::vector<double>& pts_in,
                             const std::vector<double>& pts_out)
{
  for (unsigned int i=0; i < m_inputs.size(); ++i)
    if (isSame(pts_in, pts_out, m_inputs[i]))
      return i;

  // else create new basis
  m_inputs.emplace_back(pts_in, pts_out);
  m_vals.emplace_back(boost::extents[pts_out.size()][pts_in.size()]);
  m_derivs.emplace_back(boost::extents[pts_out.size()][pts_in.size()]);
  createBasis(pts_in, pts_out, m_vals.back(), m_derivs.back());
  return m_vals.size() - 1;
}

void LagrangeMemoizer::createBasis(const std::vector<double>& pts_in,
                                   const std::vector<double>& pts_out,
                                   ArrayType<Real, 2>& vals,
                                   ArrayType<Real, 2>& derivs)
{
  LagrangeBasis basis(pts_in);
  for (Index i=0; i < pts_out.size(); ++i)
    for (Index j=0; j < pts_in.size(); ++j)
    {
      vals[i][j] = basis.evalPoly(j, pts_out[i]);
      derivs[i][j] = basis.evalPolyDeriv(j, pts_out[i]);
    }
}

bool LagrangeMemoizer::isSame(const std::vector<double>& pts_in,
                              const std::vector<double>& pts_out,
                              const InputOutput& io)
{
  if (pts_in.size() != io.first.size())
    return false;

  if (pts_out.size() != io.second.size())
    return false;

  for (unsigned int i=0; i < pts_in.size(); ++i)
    if ( std::abs(pts_in[i] - io.first[i]) > m_eps)
      return false;

  for (unsigned int i=0; i < pts_out.size(); ++i)
    if ( std::abs(pts_out[i] - io.second[i]) > m_eps)
      return false;

  return true;
}

LagrangeMemoizer lagrange_memoizer = LagrangeMemoizer();  // definition
