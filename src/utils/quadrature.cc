#include "utils/quadrature.h"

Quadrature getGaussianQuadrature(const int exactness)
{
  // exact for 2 * npts - 1,
  std::pair<Real, Real> domain = {-1, 1};
  switch (exactness)
  {
    case 1:
    {
      std::vector<Real> pts     = {0};
      std::vector<Real> weights = {2};
      int exactness = 1;
      return Quadrature(pts, weights, exactness, domain);
    }
    case 3:
    {
      std::vector<Real> pts = {-1.0/std::sqrt(3), 1.0/std::sqrt(3)};
      std::vector<Real> weights = {1, 1};
      int exactness = 3;
      return Quadrature(pts, weights, exactness, domain);
    }

    case 5:
    {
      std::vector<Real> pts     = {-std::sqrt(3.0/5.0), 0, std::sqrt(3.0/5.0)};
      std::vector<Real> weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
      int exactness = 5;
      return Quadrature(pts, weights, exactness, domain);
    }

    case 7:
    {
      Real val_tmp = (2.0/7.0) * std::sqrt(6.0/5.0);
      std::vector<Real> pts = {-std::sqrt(3.0/7.0 + val_tmp), -std::sqrt(3.0/7.0 - val_tmp),
                                std::sqrt(3.0/7.0 - val_tmp),  std::sqrt(3.0/7.0 + val_tmp)};
      std::vector<Real> weights = {(18.0 - std::sqrt(30.0))/36.0, (18.0 + std::sqrt(30.0))/36.0,
                                   (18.0 + std::sqrt(30.0))/36.0, (18.0 - std::sqrt(30.0))/36.0};
      int exactness = 7;
      return Quadrature(pts, weights, exactness, domain);
    }

    default:
      throw std::invalid_argument(std::string("no Gaussian quadrature rule for degree " + std::to_string(exactness)));

  }
}
