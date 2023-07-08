#ifndef LINEAR_SYSTEM_AUGMENTED_ASSEMBLER_H
#define LINEAR_SYSTEM_AUGMENTED_ASSEMBLER_H

#include "ProjectDefs.h"
#include "linear_system/large_matrix.h"
#include "mesh/mesh.h"
#include "mpi_utils.h"
#include <mpi.h>

namespace linear_system {

// For assembling values into a sparse matrix created with SparsityPatternAugmented
class AugmentedAssembler
{
  public:
    AugmentedAssembler(DiscPtr disc, LargeMatrixPtr mat, int num_augmented) :
      m_comm(MPI_COMM_WORLD),
      m_am_i_last_rank(commRank(m_comm) == (commSize(m_comm)-1)),
      m_num_augmented(num_augmented),
      m_mat(mat),
      m_dofs_rows(num_augmented),
      m_vals_rows(num_augmented),
      m_dofs_rows_reqs(num_augmented),
      m_vals_rows_reqs(num_augmented),
      //m_recv_dofs_rows(m_am_i_last_rank ? commSize(m_comm) : 0),
      //m_recv_vals_rows(m_am_i_last_rank ? commSize(m_comm) : 0),
      //m_recv_dofs_rows_reqs(m_am_i_last_rank ? commSize(m_comm) : 0),
      //m_recv_vals_rows_reqs(m_am_i_last_rank ? commSize(m_comm) : 0)
      m_recv_dofs_rows(commSize(m_comm)),
      m_recv_vals_rows(commSize(m_comm)),
      m_recv_dofs_rows_reqs(commSize(m_comm)),
      m_recv_vals_rows_reqs(commSize(m_comm))
    {
      disc->getMesh()->getLocalToGlobalDofs(m_local_dof_to_global); 

      //if (m_am_i_last_rank)
      {
        for (int i=0; i < commSize(m_comm); ++i)
        {
          m_recv_dofs_rows[i].resize(num_augmented);
          m_recv_vals_rows[i].resize(num_augmented);
          m_recv_dofs_rows_reqs[i].resize(num_augmented, MPI_REQUEST_NULL);
          m_recv_vals_rows_reqs[i].resize(num_augmented, MPI_REQUEST_NULL);
        }

        if (m_am_i_last_rank)
          m_recv_counts_rows.resize(commSize(m_comm)*num_augmented);
      }

      DofInt max_mesh_dof_local = *(std::max_element(m_local_dof_to_global.begin(), m_local_dof_to_global.end()));
      MPI_Allreduce(&max_mesh_dof_local, &m_num_mesh_dofs, 1, DofInt_MPI_DATATYPE, MPI_MAX, MPI_COMM_WORLD);
      m_num_mesh_dofs++;
    }

    bool amILastRank() const { return m_am_i_last_rank; }

    void setAlpha(Real alpha) { m_alpha = alpha; }

    Real getAlpha() const { return m_alpha; }

    void assembleValuesColumn(const std::vector<DofInt>& row_local_dofs, int augmented_column, const std::vector<Real>& vals)
    {
      assert(augmented_column >= 0 && augmented_column < m_num_augmented);
      assert(vals.size() == row_local_dofs.size());

      //TODO: dynamic memory allocation
      std::vector<DofInt> col_dofs = {m_num_mesh_dofs + augmented_column};
      std::vector<DofInt> row_dofs = row_local_dofs;
      ArrayType<Real, 2> vals_mat(boost::extents[row_dofs.size()][col_dofs.size()]);

      for (auto& val : row_dofs)
        val = m_local_dof_to_global[val];

      for (size_t i=0; i < row_dofs.size(); ++i)
        vals_mat[i][0] = m_alpha * vals[i];

      m_mat->assembleValues(row_dofs, col_dofs, vals_mat);
    }

    void assembleValuesRow(int augmented_row, const std::vector<DofInt>& column_local_dofs, const std::vector<Real>& vals)
    {
      assert(augmented_row >= 0 && augmented_row < m_num_augmented);
      assert(vals.size() == column_local_dofs.size());

      for (size_t i=0; i < vals.size(); ++i) 
      {
        m_dofs_rows[augmented_row].push_back(m_local_dof_to_global[column_local_dofs[i]]);
        m_vals_rows[augmented_row].push_back(m_alpha * vals[i]);
      }
    }

    void assembleAugmentedValuesDiag(const std::vector<DofInt>& augmented_rows, const std::vector<DofInt>& augmented_columns,
                                     ArrayType<Real, 2>& vals)
    {
      if (!m_am_i_last_rank)
        throw std::runtime_error("can only assemble diagonal block of augmented values on the last process");

      std::vector<DofInt> row_dofs = augmented_rows;
      std::vector<DofInt> col_dofs = augmented_columns;

      for (auto& dof : row_dofs)
        dof += m_num_mesh_dofs;

      for (auto& dof : col_dofs)
        dof += m_num_mesh_dofs;

      for (size_t i=0; i < vals.shape()[0]; ++i)
        for (int j=0; j < vals.shape()[0]; ++j)
          vals[i][j] *= m_alpha;
        
      m_mat->assembleValues(row_dofs, col_dofs, vals);      
    }

    void startAssembly()
    {
      std::cout << "\nstarting assembling" << std::endl;
      std::cout << "myrank = " << commRank(m_comm) << " / " << commSize(m_comm) << std::endl;
      getRecvCounts();
      startSends();
      startRecvs();
    }

    void finishAssembly()
    {
      std::cout << "\nfinishing assembling" << std::endl;

      finishRecvs();
      finishSends();
      for (auto& buf : m_dofs_rows)
        buf.resize(0);

      for (auto& buf : m_vals_rows)
        buf.resize(0);
    }

  private:
    void getRecvCounts()
    {
      int root = commSize(m_comm) - 1;

      std::vector<DofInt> send_counts_rows(m_num_augmented);
      for (int i=0; i < m_num_augmented; ++i)
        send_counts_rows[i] = m_vals_rows[i].size();

      if (m_am_i_last_rank)
        m_recv_counts_rows.resize(m_num_augmented*commSize(m_comm));

      MPI_Gather(send_counts_rows.data(), m_num_augmented, DofInt_MPI_DATATYPE, m_recv_counts_rows.data(),
                 m_num_augmented, MPI_INT, root, m_comm);
    }

    void startSends()
    {
      int tag = m_start_tag;
      int root = commSize(m_comm) - 1;
      
      for (int i=0; i < m_num_augmented; ++i)
      {
        MPI_Isend(m_dofs_rows[i].data(), m_dofs_rows[i].size(), DofInt_MPI_DATATYPE, root, tag++, m_comm, &(m_dofs_rows_reqs[i]));
        MPI_Isend(m_vals_rows[i].data(), m_vals_rows[i].size(), MPI_DOUBLE, root, tag++, m_comm, &(m_vals_rows_reqs[i]));
      }
    }

    void startRecvs()
    {
      if (!m_am_i_last_rank)
        return;
        
      int comm_size = commSize(m_comm);

      for (int rank=0; rank < comm_size; ++rank)
      {
        int tag = m_start_tag;
        std::cout << "rank = " << rank << std::endl;

        for (int i=0; i < m_num_augmented; ++i)
        {
          std::cout << "  i = " << i << ", number of values = " << m_recv_dofs_rows[rank][i].size() << std::endl;
          m_recv_dofs_rows[rank][i].resize(m_recv_counts_rows[m_num_augmented * rank + i]);
          m_recv_vals_rows[rank][i].resize(m_recv_counts_rows[m_num_augmented * rank + i]);


          MPI_Irecv(m_recv_dofs_rows[rank][i].data(), m_recv_dofs_rows[rank][i].size(), DofInt_MPI_DATATYPE, rank,
                    tag++, m_comm, &(m_recv_dofs_rows_reqs[rank][i]));
          std::cout << "  recv_req = " << m_recv_dofs_rows_reqs[rank][i] << std::endl;
          MPI_Irecv(m_recv_vals_rows[rank][i].data(), m_recv_vals_rows[rank][i].size(), MPI_DOUBLE, rank,
                    tag++, m_comm, &(m_recv_vals_rows_reqs[rank][i]));
        }
      }
    }

    void finishRecvs()
    {
      if (!m_am_i_last_rank)
        return;

      int comm_size = commSize(m_comm);

      for (int rank=0; rank < comm_size; ++rank)
      {
        MPI_Waitall(m_num_augmented, m_recv_dofs_rows_reqs[rank].data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(m_num_augmented, m_recv_vals_rows_reqs[rank].data(), MPI_STATUSES_IGNORE);
        std::vector<DofInt> row_dofs(1);


        for (int i=0; i < m_num_augmented; ++i)
        {
          row_dofs[0] = i + m_num_mesh_dofs;

          int nvals = m_recv_vals_rows[rank][i].size();
          ArrayType<Real, 2> vals(boost::extents[1][nvals]);
          for (int j=0; j < nvals; ++j)
            vals[0][j] = m_recv_vals_rows[rank][i][j];

          m_mat->assembleValues(row_dofs, m_recv_dofs_rows[rank][i], vals);
        }
      }
    }

    void finishSends()
    {
      MPI_Waitall(m_num_augmented, m_dofs_rows_reqs.data(), MPI_STATUSES_IGNORE);
      MPI_Waitall(m_num_augmented, m_vals_rows_reqs.data(), MPI_STATUSES_IGNORE);
    }

    MPI_Comm m_comm;
    bool m_am_i_last_rank;
    int m_num_augmented;
    DofInt m_num_mesh_dofs;
    LargeMatrixPtr m_mat;
    std::vector<DofInt> m_local_dof_to_global;
    Real m_alpha = 1.0;

    std::vector<std::vector<DofInt>> m_dofs_rows;
    std::vector<std::vector<Real>> m_vals_rows;
    std::vector<MPI_Request> m_dofs_rows_reqs;
    std::vector<MPI_Request> m_vals_rows_reqs;

    std::vector<int> m_recv_counts_rows;

    std::vector<std::vector<std::vector<DofInt>>> m_recv_dofs_rows;
    std::vector<std::vector<std::vector<Real>>>   m_recv_vals_rows;

    std::vector<std::vector<MPI_Request>> m_recv_dofs_rows_reqs;
    std::vector<std::vector<MPI_Request>> m_recv_vals_rows_reqs;

    const int m_start_tag = 1000;
};

using AugmentedAssemblerPtr = std::shared_ptr<AugmentedAssembler>;


}  // namespace

#endif