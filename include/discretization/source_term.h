#ifndef SOURCE_TERM_H
#define SOURCE_TERM_H

#include "discretization/volume_discretization.h"
//#include "physics/heat/bc_defs.h"


class SourceTerm
{
  public:
    SourceTerm (VolDiscPtr vol) : m_vol(vol) {}

    virtual ~SourceTerm() {}

    VolDiscPtr getVolDisc() const { return m_vol; }

    // gets source term values at quadrature nodes
    virtual void getValues(Index el, const Real t, Real* vals) = 0;

  protected:

    // gets coordinates of volume nodes in flat arrya
    void getQuadNodeCoords(Index el, ArrayType<Real, 2>& quad_coords)
    {
      for (int d=0; d < 3; ++d)
      {
        auto coords_d = m_vol->vol_group.coords[boost::indices[el][range()][d]];
        auto quad_coords_d = quad_coords[boost::indices[range()][d]];
        m_vol->interp_cq_flat_to_flat.interpolateVals(coords_d, quad_coords_d);
      }
    }

  private:
    VolDiscPtr m_vol;
};

using SourceTermPtr = std::shared_ptr<SourceTerm>;


// T must be callable with signature f(x, y, z, t)
template <typename T>
class SourceTermMMS : public SourceTerm
{
  public:
    SourceTermMMS(VolDiscPtr vol, T func) : 
      SourceTerm(vol),
      m_func(func),
      m_coords(boost::extents[vol->getNumQuadPtsPerElement()][3])
    {}

    void getValues(Index el, const Real t, Real* vals) override
    {
      getQuadNodeCoords(el, m_coords);
      for (unsigned int i=0; i < m_coords.shape()[0]; ++i)
        vals[i] = m_func(m_coords[i][0], m_coords[i][1], m_coords[i][2], t);
    }


  private:
    T m_func;
    ArrayType<Real, 2> m_coords;
};


template <typename T>
SourceTermPtr makeSourcetermMMS(VolDiscPtr vol, T func)
{
  return std::make_shared<SourceTermMMS<T>>(vol, func);
}




#endif