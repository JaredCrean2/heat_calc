#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <memory>
#include <vector>

#include "ProjectDefs.h"
#include "discretization/discretization.h"
#include "discretization/dof_numbering.h"
#include "large_matrix.h"

namespace linear_system {

class Assembler
{
  public:
    explicit Assembler(DiscPtr disc, LargeMatrixPtr mat) :
      m_disc(disc),
      m_dof_nums(disc->getDofNumbering()),
      m_alpha(1),
      m_vol_dofs(disc->getVolDisc(0)->getNumSolPtsPerElement()),
      m_face_dofs(disc->getSurfDisc(0)->getNumSolPtsPerFace()),
      m_matrix(mat)
    {}

    virtual ~Assembler() {}

    void setAlpha(Real alpha) { m_alpha = alpha; }

    Real getAlpha() const { return m_alpha; }

    int getNumDofs() const { return m_dof_nums->getNumDofs(); }

    // jac is df_i/du_j, where i, j in [0, numSolPtsPerElement)
    // assembla alpha * jac into the matrix
    void assembleVolume(int vol_disc_idx, int elnum, ArrayType<Real, 2>& jac);

    // jac is df_i/du_j, where i, j in [0, numSolPtsPerFace]
    void assembleFace(int surf_disc_idx, int facenum, ArrayType<Real, 2>& jac);

    void zeroMatrix() { m_matrix->zeroMatrix(); }

  private:
    DiscPtr m_disc;
    DofNumberingPtr m_dof_nums;
    Real m_alpha;
    std::vector<DofInt> m_vol_dofs;
    std::vector<DofInt> m_face_dofs;
    LargeMatrixPtr m_matrix;
};

using AssemblerPtr = std::shared_ptr<Assembler>;

} // namespace

#endif
