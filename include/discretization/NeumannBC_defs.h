#ifndef NEUMANNBC_DEFS_H
#define NEUMANNBC_DEFS_H

#include "NeumannBC.h"

template <typename T>
class NeumannBCMMS : public NeumannBC
{
  public:
    // func must be std::array<Real, 3> func(Real x, Real y, Real z, Real t)
    NeumannBCMMS(SurfDiscPtr surf, T func) :
      NeumannBC(surf, false), 
      m_func(func),
      m_coords(boost::extents[surf->getNumQuadPtsPerFace()][3])
  {}


    void getValue(const Index face, const Real t, const Real* sol_vals, Real* vals) override
    {
      auto surf_disc = getSurfDisc();
      getQuadNodeCoords(face, m_coords);

      for (int i=0; i < surf_disc->getNumQuadPtsPerFace(); ++i)
      {
        auto vals_i = m_func(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
        for (int d=0; d < 3; ++d)
          vals[i + d * surf_disc->getNumQuadPtsPerFace()] = vals_i[d];
      }
    }


  private:
      T m_func;
      ArrayType<Real, 2> m_coords;
};

template <typename T>
NeumannBCPtr makeNeumannBCMMS(SurfDiscPtr surf, T func)
{
  return std::make_shared<NeumannBCMMS<T>>(surf, func);
}


#endif