#include "linear_system/sparsity_pattern_augmented.h"
#include <iostream>
#include "mpi_utils.h"

namespace linear_system {

std::ostream& operator<<(std::ostream& os, const std::vector<PetscInt>& vec)
{
  os << "{ ";
  for (auto& val : vec)
    os << val << ", ";
  os << "}";

  return os;
}


SparsityPatternAugmented::SparsityPatternAugmented(std::shared_ptr<SparsityPattern> base_pattern, int num_augmented_rows, MPI_Comm comm) :
  m_base_pattern(base_pattern),
  m_num_augmented_rows(num_augmented_rows),
  m_comm(comm),
  m_am_i_last_rank(commRank(comm) == (commSize(comm) - 1))
{
  PetscInt num_base_dofs = base_pattern->getNumOwnedDofs();
  MPI_Allreduce(&num_base_dofs, &m_num_global_dofs_base, 1, MPIU_INT, MPI_SUM, m_comm);
}

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
    if (m_am_i_last_rank)
    {
      int num_owned_dofs_base = m_onproc_dofs.size();
      for (auto& val : m_onproc_dofs)
        val += m_num_augmented_rows;

      for (int i=0; i < m_num_augmented_rows; ++i)
        m_onproc_dofs.push_back(num_owned_dofs_base + m_num_augmented_rows);
    }
  }

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
    PetscInt num_owned_dofs = m_am_i_last_rank ? 0 : getNumOwnedDofs();
    PetscInt num_off_proc_dofs_on_last_proc = 0;
    MPI_Reduce(&num_owned_dofs, &num_off_proc_dofs_on_last_proc, 1, MPIU_INT, MPI_SUM, commSize(m_comm) - 1, m_comm);

    m_remote_dofs = m_base_pattern->getOffProcCounts();
    if (m_am_i_last_rank)
    {
      for (int i=0; i < m_num_augmented_rows; ++i)
        m_remote_dofs.push_back(num_off_proc_dofs_on_last_proc);
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
  if (!m_am_i_last_rank)
  {
    m_ghost_dofs_to_global = m_base_pattern->getGhostGlobalIndices();

    PetscInt dof = m_num_global_dofs_base;
    for (int i=0; i < m_num_augmented_rows; ++i)
      m_ghost_dofs_to_global.push_back(dof++);

    return m_ghost_dofs_to_global;
  } else
  {
    return m_base_pattern->getGhostGlobalIndices();
  }
}

const std::vector<PetscInt>& SparsityPatternAugmented::getGhostLocalIndices()
{
  if (m_ghost_dofs_to_local.size() == 0)
  {
    if (m_am_i_last_rank)
    {
      //TODO: is this right?  If the augmented dofs come after the local ones,
      //      this doesn't seem right
      m_ghost_dofs_to_local = m_base_pattern->getGhostLocalIndices();
      for (auto& val : m_ghost_dofs_to_local)
        val += m_num_augmented_rows;
    } else
    {
      m_ghost_dofs_to_local = m_base_pattern->getGhostLocalIndices();
      PetscInt dof = m_base_pattern->getNumLocalDofs();
      for (int i=0; i < m_num_augmented_rows; ++i)
        m_ghost_dofs_to_local.push_back(dof++);
    }
  }

  return m_ghost_dofs_to_local;

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

const std::vector<PetscInt>& SparsityPatternAugmented::getLocalToGlobalDofs()
{
  if (m_local_dof_to_global.size() == 0)
  {
    m_local_dof_to_global = m_base_pattern->getLocalToGlobalDofs();

    PetscInt dof = m_num_global_dofs_base;
    for (int i=0; i < m_num_augmented_rows; ++i)
      m_local_dof_to_global.push_back(dof++);    
  }


  return m_local_dof_to_global;
}


}  // namespace