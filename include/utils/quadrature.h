#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "ProjectDefs.h"
#include <stdexcept>
#include <vector>
#include <cmath>

class Quadrature
{
  public:
    Quadrature(const std::vector<Real>& pts,
               const std::vector<Real>& weights,
               const int exactness,
               const std::pair<Real, Real>& domain) :
      m_pts(pts),
      m_weights(weights),
      m_exactness(exactness),
      m_domain(domain)
    {}

    Real getPoint(const int i) const { return m_pts[i]; }

    Real getWeight(const int i) const { return m_weights[i]; }

    int getExactness() const { return m_exactness; }

    int getNumPoints() const { return m_pts.size(); }

    std::pair<Real, Real> getDomain() const { return m_domain; }

  private:
    std::vector<Real> m_pts;
    std::vector<Real> m_weights;
    int m_exactness;  // degree of polynomials this quadrature is exact for
    std::pair<Real, Real> m_domain;
};


Quadrature getGaussianQuadrature(const int exactness);

#endif
