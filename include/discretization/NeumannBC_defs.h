#ifndef NEUMANNBC_DEFS_H
#define NEUMANNBC_DEFS_H

#include "NeumannBC.h"

template <typename T>
class NeumannBCMMS : public NeumannBC
{
  public:
    NeumannBCMMS(SurfDiscPtr surf, T func) :
      NeumannBC(surf), 
      m_func(func),
      m_coords(boost::extents[surf->getNumQuadPtsPerFace()][3])
  {}


    void getValue(const Index face, const Real t, const Real* sol_vals, Real* vals) override
    {
      auto surf_disc = getSurfDisc();
      getQuadNodeCoords(face, m_coords);

      for (int i=0; i < surf_disc->getNumQuadPtsPerFace(); ++i)
        vals[i] = m_func(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
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