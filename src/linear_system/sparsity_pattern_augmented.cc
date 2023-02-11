#include "linear_system/sparsity_pattern_augmented.h"
#include <iostream>

namespace linear_system {

std::ostream& operator<<(std::ostream& os, const std::vector<PetscInt>& vec)
{
  os << "{ ";
  for (auto& val : vec)
    os << val << ", ";
  os << "}";

  return os;
}


SparsityPatternAugmented::SparsityPatternAugmented(std::shared_ptr<SparsityPattern> base_pattern, int num_augmented_rows, bool am_i_last_rank) :
  m_base_pattern(base_pattern),
  m_num_augmented_rows(num_augmented_rows),
  m_am_i_last_rank(am_i_last_rank)
{}

PetscInt SparsityPatternAugmented::getNumOwnedDofs() const
{
  return m_base_pattern->getNumOwnedDofs() + (m_am_i_last_rank ? m_num_augmented_rows : 0);
}

// returns vector giving number of local dofs connected to each dof
const std::vector<PetscInt>& SparsityPatternAugmented::getDiagonalCounts()
{
  if (m_onproc_dofs.size() == 0)
  {
    m_onproc_dofs = m_base_pattern->getDiagonalCounts();
    std::cout << "base onproc dofs = " << m_onproc_dofs << std::endl;
    if (m_am_i_last_rank)
    {
      std::cout << "I am the last rank" << std::endl;
      int num_owned_dofs_base = m_onproc_dofs.size();
      for (auto& val : m_onproc_dofs)
        val += m_num_augmented_rows;

      for (int i=0; i < m_num_augmented_rows; ++i)
        m_onproc_dofs.push_back(num_owned_dofs_base + m_num_augmented_rows);
    }
  }

  std::cout << "augmented dofs = " << m_onproc_dofs << std::endl;
  return m_onproc_dofs;
}

// similar to the above, but returns the count for only the matrix diagonal + upper triangle
const std::vector<PetscInt>& SparsityPatternAugmented::getDiagonalCountsSym()
{
  throw std::runtime_error("symmetric matrix not supported with augmented rows");
}

// returns vector giving number of remote dofs connected to each dof
const std::vector<PetscInt>& SparsityPatternAugmented::getOffProcCounts()
{
  if (m_remote_dofs.size() == 0)
  {
    m_remote_dofs = m_base_pattern->getOffProcCounts();
    if (m_am_i_last_rank)
    {
      for (int i=0; i < m_num_augmented_rows; ++i)
        m_remote_dofs.push_back(0);
    } else    
    {
      for (auto& val : m_remote_dofs)
        val += m_num_augmented_rows;
    }
  }

  return m_remote_dofs;
}
    
// similar to the above, but returns the count for only the upper triangle
const std::vector<PetscInt>& SparsityPatternAugmented::getOffProcCountsSym()
{
  throw std::runtime_error("symmetric matrix not supported with augmented rows");
}

const std::vector<PetscInt>& SparsityPatternAugmented::getGhostGlobalIndices()
{
  return m_base_pattern->getGhostGlobalIndices();
}

const std::vector<PetscInt>& SparsityPatternAugmented::getGhostLocalIndices()
{
  return m_base_pattern->getGhostLocalIndices();
}

const std::vector<PetscInt>& SparsityPatternAugmented::getOwnedToLocalInfo()
{
  if (m_am_i_last_rank)
  {
    if (m_owned_dof_to_local.size() == 0)
    {
      m_owned_dof_to_local = m_base_pattern->getOwnedToLocalInfo();
      int num_local_dofs_base = m_owned_dof_to_local.size() + getGhostLocalIndices().size();
      for (int i=0; i < m_num_augmented_rows; ++i)
        m_owned_dof_to_local.push_back(num_local_dofs_base + i);
    }

    return m_owned_dof_to_local;

  } else
  {
    return m_base_pattern->getOwnedToLocalInfo();
  }
}

}  // namespace