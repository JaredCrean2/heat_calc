#ifndef DISCVECTOR_H
#define DISCVECTOR_H

#include "ProjectDefs.h"
#include "discretization/discretization.h"
#include "discretization/dof_numbering.h"

#include <stdexcept>

class DiscVector
{
  public:
    explicit DiscVector(DiscPtr disc);

    void markVectorModified();

    void markArrayModified();

    // overwrites array with vector.  Does not overwrite all values because
    // Dirichlet values are not present in the vector
    void syncVectorToArray();

    // sums array entries into vector.
    //   zero_vec: if true (default) zeros the vector first
    void syncArrayToVector(const bool zero_vec = true);

    ArrayType<Real, 1>& getVector();

    ArrayType<Real, 2>& getArray(const Index idx);

    // sets both vector and array (including dirchlet values) to a
    // constant value
    void set(const Real val);

    // sets both the vector and array (including dirichlet values) to
    // a given function f(x,y,z)
    // The function is evaluated at the solution nodes
    template <typename T>
    void setFunc(T func);

    DiscPtr getDisc() const { return m_disc;}

  private:
    DiscPtr m_disc;
    ArrayType<Real, 1> m_vec;    // length numdofs
    std::vector<ArrayType<Real, 2>> m_array;  // vector length num volume
                                              // groups,
                                              // inner array nelem per group x numSolNodesPerElement
    bool m_is_vec_modified   = false;
    bool m_is_array_modified = false;
};


template <typename T>
void DiscVector::setFunc(T func)
{
  auto dof_numbering = m_disc->getDofNumbering();
  for (int i=0; i < m_disc->getNumVolDiscs(); ++i)
  {
    VolDiscPtr disc = m_disc->getVolDisc(i);
    const auto& dof_nums = dof_numbering->getDofs(i);
    auto& array_i   = m_array[i];
    //TODO: use one stored in disc
    LagrangeEvaluatorTPFlatToTPFlat mapper(
                         disc->vol_group.getTPMapperCoord().getXi(),
                         disc->vol_group.getTPMapperSol().getXi(),
                         disc->vol_group.getTPMapperCoord().getNodemap(),
                         disc->vol_group.getTPMapperSol().getNodemap());
    ArrayType<Real, 2> coords_sol(boost::extents[disc->getNumSolPtsPerElement()][3]);
    for (int el=0; el < disc->getNumElems(); ++el)
    {
      for (int d=0; d < 3; ++d)
      {
        auto coords_el = disc->vol_group.coords[boost::indices[el][range()][d]];
        auto coords_d = coords_sol[boost::indices[range()][d]];
        mapper.interpolateVals(coords_el, coords_d);
      }

      for (int node=0; node < disc->getNumSolPtsPerElement(); ++node)
      {
        auto val = func(coords_sol[node][0], coords_sol[node][1],
                        coords_sol[node][2]);
        auto dof = dof_nums[el][node];
        array_i[el][node] = val;
        if (dof_numbering->isDofActive(dof))
          m_vec[dof] = val;
      }
    }
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
};

using DiscVectorPtr = std::shared_ptr<DiscVector>;

DiscVectorPtr makeDiscVector(DiscPtr disc);

#endif
