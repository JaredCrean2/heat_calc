#ifndef DIRICHLETBC_DEFS_H
#define DIRICHLETBC_DEFS_H

#include "discretization/DirichletBC.h"

template <typename T>
class DirichletBCMMS : public DirichletBC
{
  public:
    DirichletBCMMS(SurfDiscPtr surf, T func) : m_func(func) {}


    void getValue(const Index face, const Real t, Real* vals) const override
    {
      auto surf_disc = getSurfDisc();
      ArrayType<Real, 2> coords(boost::extents[surf_disc->getNumQuadPointPtsPerFace()][3]); //TODO: store this
      getNodeCoords(face, coords);

      for (int i=0; i < getSurfDisc()->getNumCoordPtsPerFace(); ++i)
        vals[i] = m_func(coords[i][0], coords[i][1], coords[i][2], t);
    }

  private:
      T m_func;
};

template <typename T>
DirichletBCPtr makeDirichletBCMMS(SurfDiscPtr surf, T func)
{
  return std::make_shared<DirichletBCMMS<T>>(surf, func);
}

#endif
