#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <memory>
#include <vector>

//TODO: could forward declare stuff from discretization and dof_numbering.h
#include "ProjectDefs.h"
#include "discretization/discretization.h"
#include "discretization/dof_numbering.h"
#include "large_matrix.h"

namespace linear_system {

class Assembler
{
  public:
    explicit Assembler(DiscPtr disc, LargeMatrixPtr mat);

    virtual ~Assembler() {}

    void setAlpha(Real alpha) { m_alpha = alpha; }

    Real getAlpha() const { return m_alpha; }

    int getNumLocalDofs() const { return m_dof_nums->getNumLocalDofs(); }

    int getNumOwnedDofs() const { return m_dof_nums->getNumOwnedDofs(); }

    // jac is df_i/du_j, where i, j in [0, numSolPtsPerElement)
    // assembla alpha * jac into the matrix
    void assembleVolume(int vol_disc_idx, int elnum, ArrayType<Real, 2>& jac);

    // jac is df_i/du_j, where i, j in [0, numSolPtsPerFace]
    void assembleFace(int surf_disc_idx, int facenum, ArrayType<Real, 2>& jac);

    void zeroMatrix() { m_matrix->zeroMatrix(); }

    LargeMatrixPtr getMatrix() { return m_matrix; }

  private:
    DiscPtr m_disc;
    DofNumberingPtr m_dof_nums;
    std::vector<DofInt> m_local_dof_to_global;
    Real m_alpha;
    std::vector<DofInt> m_vol_dofs;
    std::vector<DofInt> m_face_dofs;
    LargeMatrixPtr m_matrix;
};

using AssemblerPtr = std::shared_ptr<Assembler>;

} // namespace

#endif
