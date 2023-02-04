#ifndef TIMESTEP_CONTROLLER_H
#define TIMESTEP_CONTROLLER_H

#include "ProjectDefs.h"

#include <iostream>  //TODO: DEBUGGING

namespace timesolvers {

class TimestepController 
{
  public:
    virtual ~TimestepController() {}

    virtual Real getNextTimestep(Real t) = 0;

    virtual void recordLastIteration(Real physics_rhs) = 0;
};

class TimestepControllerConstant : public TimestepController
{
  public:
    explicit TimestepControllerConstant(Real delta_t) :
      m_delta_t(delta_t)
    {}

    Real getNextTimestep(Real t) override { return m_delta_t; }

    void recordLastIteration(Real physics_rhs) override {};

  private:
    Real m_delta_t;
};


// update the timestep based on the residual norm.
// This is useful for solving problems to steady state, increasing
// the timestep as the residual converges
class TimestepControllerResidual : public TimestepController
{
  public:
    TimestepControllerResidual(Real delta_t0, Real exponent) :
      m_delta_t_prev(delta_t0),
      m_exponent(exponent)
    {}

    Real getNextTimestep(Real t) override { return m_delta_t_prev; }

    void recordLastIteration(Real physics_rhs) override
    {
      std::cout << "recording residual " << physics_rhs << std::endl;
      m_residual_prev = m_residual_curr;
      m_residual_curr = physics_rhs;

      if (m_residual_prev != -1)
      {
        Real ratio = std::pow(m_residual_curr/m_residual_prev, m_exponent);
        m_delta_t_prev /= ratio;
        std::cout << "new delta t = " << m_delta_t_prev << std::endl;
      }
    };


  private:
    Real m_delta_t_prev;
    Real m_exponent;
    Real m_residual_prev = -1;
    Real m_residual_curr = -1;


};

}  // namespace

#endif