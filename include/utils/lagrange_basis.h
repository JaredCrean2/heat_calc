#ifndef LAGRANGE_BASIS_H
#define LAGRANGE_BASIS_H

#include "ProjectDefs.h"
#include <iostream>

#include <utility>

// class for evaluating Lagrange polynomials and their derivatives
class LagrangeBasis
{
  public: 
    using Index = std::vector<double>::size_type;

    LagrangeBasis(const std::vector<double>& pts);
    // evaluate the pth polynomial at point x
    double evalPoly(Index p, double x);

    double evalPolyDeriv(Index p, double x);

    Index getNumPoints() const { return m_pts.size();}

  private:
    std::vector<double> m_pts;
    std::vector<double> m_denoms;

};


class LagrangeMemoizer
{
  public:
    using SType   = ArrayType<Real, 2>;
    using RetType = const SType&;
    using RetTypeP = SType*;

    RetType getValues(const std::vector<double>& pts_in,
                                        const std::vector<double>& pts_out)
    {
      return m_vals.at(getIdx(pts_in, pts_out));
    }

    RetTypeP getValuesP(const std::vector<double>& pts_in,
                                        const std::vector<double>& pts_out)
    {
      return &(m_vals.at(getIdx(pts_in, pts_out)));
    }

    RetType getDerivs(const std::vector<double>& pts_in,
                                        const std::vector<double>& pts_out)
    {
      return m_derivs.at(getIdx(pts_in, pts_out));
    }

    RetTypeP getDerivsP(const std::vector<double>& pts_in,
                                        const std::vector<double>& pts_out)
    {
      return &(m_derivs.at(getIdx(pts_in, pts_out)));
    }

  private:

    using InputOutput = std::pair<std::vector<double>, std::vector<double>>;
    using Index = LagrangeBasis::Index;

    int getIdx(const std::vector<double>& pts_in,
               const std::vector<double>& pts_out);

    void createBasis(const std::vector<double>& pts_in,
                     const std::vector<double>& pts_out,
                     ArrayType<Real, 2>& vals,
                     ArrayType<Real, 2>& derivs);

    bool isSame(const std::vector<double>& pts_in,
                const std::vector<double>& pts_out,
                const InputOutput& io);

    double m_eps = 1e-13;
    std::vector<InputOutput> m_inputs;
    std::vector<SType> m_vals;  // pts_out x pts_in
    std::vector<SType> m_derivs;  // pts_out x pts_in
};

extern LagrangeMemoizer lagrange_memoizer;
#endif
