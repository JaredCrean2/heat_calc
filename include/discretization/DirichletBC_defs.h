#ifndef DIRICHLETBC_DEFS_H
#define DIRICHLETBC_DEFS_H

#include "discretization/DirichletBC.h"

template <typename T, typename T2>
class DirichletBCMMS : public DirichletBC
{
  public:
    DirichletBCMMS(SurfDiscPtr surf, T func) :
      DirichletBC(surf, false), 
      m_func(func),
      m_coords(boost::extents[surf->getNumSolPtsPerFace()][3])
  {}

    DirichletBCMMS(SurfDiscPtr surf, T func, T2 func_dt) :
      DirichletBC(surf, true), 
      m_func(func),
      m_func_dt(func_dt),
      m_coords(boost::extents[surf->getNumSolPtsPerFace()][3])
  {}


    void getValue(const Index face, const Real t, Real* vals) override
    {
      auto surf_disc = getSurfDisc();
      getSolNodeCoords(face, m_coords);

      for (int i=0; i < surf_disc->getNumSolPtsPerFace(); ++i)
        vals[i] = m_func(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
    }

    void getValueDt(const Index face, const Real t, Real* vals) override
    {
      auto surf_disc = getSurfDisc();
      getSolNodeCoords(face, m_coords);

      for (int i=0; i < surf_disc->getNumSolPtsPerFace(); ++i)
        vals[i] = m_func_dt(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
    }

  private:
      T m_func;
      T2 m_func_dt;
      ArrayType<Real, 2> m_coords;
};

template <typename T>
DirichletBCPtr makeDirichletBCMMS(SurfDiscPtr surf, T func)
{
  return std::make_shared<DirichletBCMMS<T, impl::ErrorFuncType>>(surf, func, &impl::errorFunc);
}

template <typename T, typename T2>
DirichletBCPtr makeDirichletBCMMS(SurfDiscPtr surf, T func, T2 func_dt)
{
  return std::make_shared<DirichletBCMMS<T, T2>>(surf, func, func_dt);
}

#endif
