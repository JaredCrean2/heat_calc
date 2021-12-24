#ifndef DISCVECTOR_H
#define DISCVECTOR_H

#include "ProjectDefs.h"
#include "discretization/discretization.h"
#include "discretization/dof_numbering.h"

#include <stdexcept>
#include <functional>

class DiscVector
{
  public:
    explicit DiscVector(DiscPtr disc);

    void markVectorModified();

    void markArrayModified();

    // returns true if the most current data is in the vector
    bool isVectorCurrent() const { return !m_is_array_modified;}

    // returns true if the most current data is in the array
    bool isArrayCurrent() const { return !m_is_vec_modified;}

    // overwrites array with vector.  Does not overwrite all values because
    // Dirichlet values are not present in the vector
    void syncVectorToArray();

    // sums array entries into vector.
    //   func: a binary function z = func(x, y).  x is always the value currently
    //         in the vector (at that point in the algorithm), y is the value in
    //         the array
    //   zero_vec: if true (default) zeros the vector first
    template <typename Func=std::plus<Real>>
    void syncArrayToVector(const Func func = std::plus<Real>(), const bool zero_vec = true);

    ArrayType<Real, 1>& getVector();

    const ArrayType<Real, 1>& getVector() const;

    int getNumArrays() const { return m_array.size(); }

    ArrayType<Real, 2>& getArray(const Index idx);

    const ArrayType<Real, 2>& getArray(const Index idx) const;

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


// sums array entries into vector.
//   zero_vec: if true (default) zeros the vector first
template <typename Func>
void DiscVector::syncArrayToVector(Func func, const bool zero_vec)
{
  if (!m_is_array_modified)
  {
    std::cerr << "Warning: called syncArrayToVector when array has not been modified, vector will not be updated" << std::endl;
    return;
  }

  if (zero_vec)
    for( auto& val : m_vec)
      val = 0;

  auto dof_numbering = m_disc->getDofNumbering();
  for (int i=0; i < m_disc->getNumVolDiscs(); ++i)
  {
    VolDiscPtr vol_disc = m_disc->getVolDisc(i);
    const auto& dof_nums = dof_numbering->getDofs(i);
    auto& array_i = m_array[i];

    for (int el=0; el < vol_disc->getNumElems(); ++el)
      for (int node=0; node < vol_disc->getNumSolPtsPerElement(); ++node)
      {
        auto dof = dof_nums[el][node];
        if (dof_numbering->isDofActive(dof))
          m_vec[dof] = func(m_vec[dof], array_i[el][node]);
      }
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
}


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
