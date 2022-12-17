#ifndef UTILS_SPLINE_H
#define UTILS_SPLINE_H

#include "ProjectDefs.h"

#include <iostream>

class SingleCubicSpline
{
  public:

    void setupSpline(Real x_lower, Real x_upper, Real f_lower, Real fprime_lower, Real f_upper, Real fprime_upper)
    {
      m_x_lower = x_lower;
      m_x_upper = x_upper;
      computeCoeffs(f_lower, fprime_lower, f_upper, fprime_upper);
    }

    void setupSplineRev(Real x_lower, Real x_upper, const std::array<Real, 4>& coeffs_bar,
                        Real& f_lower_bar, Real& fprime_lower_bar, Real& f_upper_bar, Real& fprime_upper_bar)
    {
      m_x_lower = x_lower;
      m_x_upper = x_upper;
      computeCoeffsRev(coeffs_bar, f_lower_bar, fprime_lower_bar, f_upper_bar, fprime_upper_bar);
    }

    bool inRange(Real x) { return x >= m_x_lower && x <= m_x_upper; }

    Real eval(Real x)
    {
      Real t = (x - m_x_lower)/(m_x_upper - m_x_lower);

      Real term1 = m_coeffs[2] + m_coeffs[3]*t;
      Real term2 = m_coeffs[1] + t*term1;
      Real term3 = m_coeffs[0] + t*term2;

      return term3;
    }

    Real evalDot(Real x, Real x_dot, Real& result_dot)
    {
      Real t = (x - m_x_lower)/(m_x_upper - m_x_lower);
      Real t_dot = x_dot/(m_x_upper - m_x_lower);

      Real term1 = m_coeffs[2] + m_coeffs[3]*t;
      Real term1_dot = m_coeffs[3]*t_dot;

      Real term2 = m_coeffs[1] + t*term1;
      Real term2_dot = t_dot * term1 + t * term1_dot;

      Real term3 = m_coeffs[0] + t*term2;
      Real term3_dot = t_dot * term2 + t * term2_dot;

      result_dot = term3_dot;

      return term3;
    }

    Real evalRev(Real x, Real result_bar, Real& x_bar, std::array<Real, 4>& coeffs_bar)
    {
      assert(inRange(x));

      Real t = (x - m_x_lower)/(m_x_upper - m_x_lower);

      Real term1 = m_coeffs[2] + m_coeffs[3]*t;
      Real term2 = m_coeffs[1] + t*term1;
      Real term3 = m_coeffs[0] + t*term2;

      //-----------------------------------------
      for (int i=0; i < 4; ++i)
        coeffs_bar[i] = 0;

      Real term3_bar = result_bar;

      coeffs_bar[0] += term3_bar;
      Real t_bar     = term2*term3_bar;
      Real term2_bar = t*term3_bar;

      coeffs_bar[1] += term2_bar;
      t_bar += term1 * term2_bar;
      Real term1_bar =  t * term2_bar;

      coeffs_bar[2] += term1_bar;
      coeffs_bar[3] += t*term1_bar;
      t_bar += m_coeffs[3] * term1_bar;
      x_bar = t_bar/(m_x_upper - m_x_lower);


      return term3;
    }

  private:
    void computeCoeffs(Real f_lower, Real fprime_lower, Real f_upper, Real fprime_upper)
    {
      // convert df/dx to df/dt
      fprime_lower = fprime_lower * (m_x_upper - m_x_lower);
      fprime_upper = fprime_upper * (m_x_upper - m_x_lower);

      m_coeffs[0] = f_lower;
      m_coeffs[1] = fprime_lower;
      m_coeffs[2] = 3*(f_upper - f_lower)  - 2*fprime_lower - fprime_upper;
      m_coeffs[3] = fprime_upper + fprime_lower - 2*(f_upper - f_lower);
    }

    // This is a linear function, so no need for the input arguments
    void computeCoeffsRev(const std::array<Real, 4>& coeffs_bar,
                          Real& f_lower_bar, Real& fprime_lower_bar, Real& f_upper_bar, Real& fprime_upper_bar)
    {
      //Real fprime_lower2 = fprime_lower * (m_x_upper - m_x_lower);
      //Real fprime_upper2 = fprime_upper * (m_x_upper - m_x_lower);

      f_lower_bar      = 0;
      fprime_lower_bar = 0;
      f_upper_bar      = 0;
      fprime_upper_bar = 0;

      fprime_upper_bar += coeffs_bar[3];
      fprime_lower_bar += coeffs_bar[3];
      f_upper_bar      -= 2*coeffs_bar[3];
      f_lower_bar      += 2*coeffs_bar[3];

      f_upper_bar      += 3*coeffs_bar[2];
      f_lower_bar      -= 3*coeffs_bar[2];
      fprime_lower_bar -= 2*coeffs_bar[2];
      fprime_upper_bar -= coeffs_bar[2];

      fprime_lower_bar += coeffs_bar[1];

      f_lower_bar      += coeffs_bar[0];

      fprime_lower_bar = fprime_lower_bar * (m_x_upper - m_x_lower);
      fprime_upper_bar = fprime_upper_bar * (m_x_upper - m_x_lower);
    }
    

    Real m_x_lower;
    Real m_x_upper;
    std::array<Real, 4> m_coeffs;
};

#endif