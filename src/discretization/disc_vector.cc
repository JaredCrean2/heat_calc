
#include "discretization/disc_vector.h"

DiscVector::DiscVector(DiscPtr disc) : 
  m_disc(disc),
  m_vec(boost::extents[m_disc->getMesh()->getNumDofs()]),
  m_array(m_disc->getMesh()->getNumVolumeGroups())
{
  for (int i=0; i < m_disc->getMesh()->getNumVolumeGroups(); ++i)
  {
    auto& vol_group = m_disc->getMesh()->getElements(i);
    m_array[i].resize(
      boost::extents[vol_group.getNumElems()][vol_group.getNumSolPtsPerElement()]);
  }
}


void DiscVector::markVectorModified()
{ 
  if (m_is_array_modified)
    throw std::runtime_error("cannot mark vector modified when array is already modified");

  m_is_vec_modified = true;
}


void DiscVector::markArrayModified()
{ 
  if (m_is_vec_modified)
    throw std::runtime_error("cannot mark array modified when vector is already modified");

  m_is_array_modified = true;
}


// overwrites array with vector.  Does not overwrite all values because
// Dirichlet values are not present in the vector
void DiscVector::syncVectorToArray()
{
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
          array_i[el][node] = m_vec[dof];
      }
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
}


ArrayType<Real, 1>& DiscVector::getVector()
{
  if (m_is_array_modified)
    throw std::runtime_error("cannot get vector when array is already modified");

  return m_vec;
}

const ArrayType<Real, 1>& DiscVector::getVector() const
{
  if (m_is_array_modified)
    throw std::runtime_error("cannot get vector when array is already modified");

  return m_vec;
}


ArrayType<Real, 2>& DiscVector::getArray(const Index idx)
{
  if (m_is_vec_modified)
    throw std::runtime_error("cannot get array when vector is already modified");

  return m_array.at(idx);
}

const ArrayType<Real, 2>& DiscVector::getArray(const Index idx) const
{
  if (m_is_vec_modified)
    throw std::runtime_error("cannot get array when vector is already modified");

  return m_array.at(idx);
}


void DiscVector::set(const Real val)
{
  for (auto& val_v : m_vec)
    val_v = val;

  for (int i=0; i < m_disc->getNumVolDiscs(); ++i)
  {
    auto& array_i = m_array[i];
    for (unsigned int i=0; i < array_i.shape()[0]; ++i)
      for (unsigned int j=0; j < array_i.shape()[1]; ++j)
        array_i[i][j] = val;
  }

  m_is_vec_modified   = false;
  m_is_array_modified = false;
}


DiscVectorPtr makeDiscVector(DiscPtr disc)
{
  return std::make_shared<DiscVector>(disc);
}


