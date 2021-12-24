#ifndef DIRICHLETBC_DEFS_H
#define DIRICHLETBC_DEFS_H

#include "discretization/DirichletBC.h"

template <typename T>
class DirichletBCMMS : public DirichletBC
{
  public:
    DirichletBCMMS(SurfDiscPtr surf, T func) :
      DirichletBC(surf), 
      m_func(func),
      m_coords(boost::extents[surf->getNumSolPtsPerFace()][3])
  {}


    void getValue(const Index face, const Real t, Real* vals) override
    {
      auto surf_disc = getSurfDisc();
      getSolNodeCoords(face, m_coords);

      for (int i=0; i < surf_disc->getNumCoordPtsPerFace(); ++i)
        vals[i] = m_func(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
    }

  private:
      T m_func;
      ArrayType<Real, 2> m_coords;
};

template <typename T>
DirichletBCPtr makeDirichletBCMMS(SurfDiscPtr surf, T func)
{
  return std::make_shared<DirichletBCMMS<T>>(surf, func);
}

#endif
