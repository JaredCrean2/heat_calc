#ifndef PHYSICS_AUXILIARY_EQUATIONS_H
#define PHYSICS_AUXILIARY_EQUATIONS_H

#include "discretization/discretization.h"
#include "discretization/disc_vector.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_factory.h"
#include "linear_system/sparsity_pattern_dense.h"

class AuxiliaryEquationStorage;
using AuxiliaryEquationsStoragePtr = std::shared_ptr<AuxiliaryEquationStorage>;


class AuxiliaryEquationsJacobians;
using AuxiliaryEquationsJacobiansPtr = std::shared_ptr<AuxiliaryEquationsJacobians>;

// this class describes any auxiliary equations (equations other than the finite element problem)
// that the PhysicsModel needs solved.  The intent is that the equations will 
// be solved as a coupled system, for example
//  [A B][x] = [a]
//  [C D][y]   [b]
//
// where A is the Jacobian of the finite element problem, there is a single
// set of auxiliary equations with Jacobian D.  B and C are the cross-term
// Jacobians.  B is the derivative of the finite element problem wrt the
// auxiliary variables, and C is the derivative of the auxiliary equations
// wrt the finite element solution.
class AuxiliaryEquations
{
  public:
    explicit AuxiliaryEquations(DiscPtr disc) :
      m_disc(disc)
    {}

    virtual ~AuxiliaryEquations() {};

    // return the number of sets of equations, including the main finite element element problem
    virtual int getNumBlocks() const { return getNumAuxiliaryBlocks() + 1; }

    // returns the number of variables in the given block
    virtual int getBlockSize(int block) const
    {
      if (block == 0)
      {
        // Note: returning the number of local dofs really is the correct
        // thing here.  In all cases, this function returns the size a
        // vector should be on this process.
        return m_disc->getDofNumbering()->getNumLocalDofs();
      } else
        return getAuxiliaryBlockSize(block - 1);
    }

    virtual void computeRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs)
    {
      assertAlways(block != 0, "Cannot compute finite element rhs via AuxiliaryEquations");
      assertAlways(rhs.shape()[0] == getBlockSize(block), "rhs vector size is incorrect");

      computeAuxiliaryRhs(block-1, u_vec, u_aux_vec, t, rhs);
    }

    virtual void computeJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, linear_system::SimpleAssemblerPtr mat)
    {
      assertAlways(block != 0, "Cannot compute finite element Jacobian via AuxiliaryEquations");
      assertAlways(mat->getNumDofs() == getBlockSize(block), "matrix size is incorrect");
      assertAlways(mat->getNumDofs() == getBlockSize(block), "matrix size is incorrect");

      computeAuxiliaryJacobian(block-1, u_vec, u_aux_vec, t,  mat);
    }

    virtual void computeMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr mat)
    {
      assertAlways(block != 0, "Cannot compute finite element Jacobian via AuxiliaryEquations");
      assertAlways(mat->getNumDofs() == getBlockSize(block), "matrix size is incorrect");
      assertAlways(mat->getNumDofs() == getBlockSize(block), "matrix size is incorrect");   

      computeAuxiliaryMassMatrix(block-1, t,  mat);
    }

    virtual void multiplyMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
    {
      assertAlways(block != 0, "Cannot compute finite element Jacobian via AuxiliaryEquations");
      assertAlways(x.shape()[0] == getBlockSize(block), "matrix size is incorrect");
      assertAlways(b.shape()[0] == getBlockSize(block), "matrix size is incorrect"); 

      multiplyAuxiliaryMassMatrix(block-1, t, x, b);
    }

    virtual void multiplyOffDiagonal(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b)
    {
      assertAlways(x.shape()[0] == getBlockSize(jblock), "x size is incorrect");
      assertAlways(b.shape()[0] == getBlockSize(iblock), "b size is incorrect");

      if (iblock == 0)
        computeFiniteElementJacobianVectorProduct(jblock, u_vec, u_aux_vec, t, x, b);
      else
        computeAuxiliaryJacobianVectorProduct(iblock, jblock, u_vec, u_aux_vec, t, x, b);
    }

    virtual AuxiliaryEquationsJacobiansPtr getJacobians() = 0;

  protected:
    // return number of auxiliary sets of equations
    virtual int getNumAuxiliaryBlocks() const = 0;
    
    // returns number of variables in each block
    virtual int getAuxiliaryBlockSize(int block) const = 0;

    // each auxiliary block must be of the form du/dt = rhs(u, t).  This function computes the rhs
    virtual void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs) = 0;

    virtual void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    virtual void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

    // compute the diagonal Jacobian block for the given block
    virtual void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    // compute the diagonal Jacobian block for the given block
    //virtual void computeAuxiliaryJacobian(int block, Real t, linear_system::SimpleAssemblerPtr mat) = 0;

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    virtual void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    virtual void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) = 0;

  private:
    DiscPtr m_disc;
};

using AuxiliaryEquationsPtr = std::shared_ptr<AuxiliaryEquations>;


class AuxiliaryEquationStorage
{
  public:
    explicit AuxiliaryEquationStorage(AuxiliaryEquationsPtr aux_eqns) :
    m_aux_eqns(aux_eqns),
    m_vectors(0)
    {
      for (int i=0; i < aux_eqns->getNumBlocks()-1; ++i)
      {
        int size = aux_eqns->getBlockSize(i+1);
        m_vectors.emplace_back(boost::extents[size]);
      }
    }

    ArrayType<Real, 1>& getVector(int block)
    {
      assertAlways(block > 0 && block < m_aux_eqns->getNumBlocks(), "block is out of range");
      return m_vectors[block-1];
    }

    const ArrayType<Real, 1>& getVector(int block) const
    {
      assertAlways(block > 0 && block < m_aux_eqns->getNumBlocks(), "block is out of range");
      return m_vectors[block-1];
    }
    
  private:
    AuxiliaryEquationsPtr m_aux_eqns;
    std::vector<ArrayType<Real, 1>> m_vectors;
};

inline AuxiliaryEquationsStoragePtr makeAuxiliaryEquationsStorage(AuxiliaryEquationsPtr aux_eqns)
{
  return std::make_shared<AuxiliaryEquationStorage>(aux_eqns);
}


class AuxiliaryEquationsJacobians
{
  public:
    explicit AuxiliaryEquationsJacobians(int num_blocks) :
      m_num_blocks(num_blocks)
    {}

    virtual ~AuxiliaryEquationsJacobians() {}

    linear_system::LargeMatrixPtr getMatrix(int block)
    {
      assertAlways(block > 0 && block < m_num_blocks, "block is out of range");
      return getAuxiliaryBlockMatrix(block-1);
    }

    const linear_system::LargeMatrixPtr getMatrix(int block) const
    {
      assertAlways(block > 0 && block < m_num_blocks, "block is out of range");
      return getAuxiliaryBlockMatrix(block-1);
    }


  protected:
    virtual linear_system::LargeMatrixPtr getAuxiliaryBlockMatrix(int block) = 0;

    virtual const linear_system::LargeMatrixPtr getAuxiliaryBlockMatrix(int block) const = 0;

  private:
    int m_num_blocks;
};


class AuxiliaryEquationsJacobiansDense : public AuxiliaryEquationsJacobians
{
  public:
    explicit AuxiliaryEquationsJacobiansDense(AuxiliaryEquations& aux_eqns) :
      AuxiliaryEquationsJacobians(aux_eqns.getNumBlocks()),
      m_jacs(aux_eqns.getNumBlocks()-1) 
    {
      for (int i=0; i < aux_eqns.getNumBlocks()-1; ++i)
      {
        int size = aux_eqns.getBlockSize(i+1);
        auto matrix_opts = std::make_shared<linear_system::LargeMatrixOpts>();
        auto sparsity = std::make_shared<linear_system::SparsityPatternDense>(size);
        m_jacs[i] = linear_system::largeMatrixFactory(linear_system::LargeMatrixType::Dense, matrix_opts, sparsity);
      }
    }

    AuxiliaryEquationsJacobiansDense() :
      AuxiliaryEquationsJacobians(1)
    {}


  protected:
    virtual linear_system::LargeMatrixPtr getAuxiliaryBlockMatrix(int block) { return m_jacs[block]; }

    virtual const linear_system::LargeMatrixPtr getAuxiliaryBlockMatrix(int block) const { return m_jacs[block]; }

  private:
    std::vector<linear_system::LargeMatrixPtr> m_jacs;
};



class AuxiliaryEquationsNone : public AuxiliaryEquations
{
  public:
    explicit AuxiliaryEquationsNone(DiscPtr disc) :
      AuxiliaryEquations(disc),
      m_disc(disc),
      m_jacs(std::make_shared<AuxiliaryEquationsJacobiansDense>())
    {}

  protected:
    virtual int getNumAuxiliaryBlocks() const override { return 0; }
    
    // returns number of variables in each block
    virtual int getAuxiliaryBlockSize(int block) const override { return 0; }


    virtual AuxiliaryEquationsJacobiansPtr getJacobians() override
    {
      return m_jacs;
    }

    // each auxiliary block must be of the form du/dt = rhs(u, t).  This function computes the rhs
    void computeAuxiliaryRhs(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, ArrayType<Real, 1>& rhs) override
    {
      assertAlways(false, "cannot compute auxiliary rhs for AuxiliaryEquationsNone"); 
    }

    virtual void computeAuxiliaryMassMatrix(int block, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      assertAlways(false, "cannot compute mass matrix for AuxiliaryEquationsNone"); 
    }

    virtual void multiplyAuxiliaryMassMatrix(int block, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      assertAlways(false, "cannot multiply auxiliary mass matrix for AuxiliaryEquationsNone"); 
    }

    // compute the diagonal Jacobian block for the given block
    virtual void computeAuxiliaryJacobian(int block, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, linear_system::SimpleAssemblerPtr mat) override
    {
      assertAlways(false, "cannot compute auxiliary Jacobian for AuxiliaryEquationsNone"); 
    }

    // compute the Jacobian-vector product for the block the couples the finite element problem to auxiliary block jblock
    virtual void computeFiniteElementJacobianVectorProduct(int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      assertAlways(false, "cannot compute finite element Jacobian vector product for AuxiliaryEquationsNone"); 
    }

    // compute the Jacobian-vector product for the block that couples auxiliary block i to auxiliary block j
    virtual void computeAuxiliaryJacobianVectorProduct(int iblock, int jblock, DiscVectorPtr u_vec, AuxiliaryEquationsStoragePtr u_aux_vec, Real t, const ArrayType<Real, 1>& x, ArrayType<Real, 1>& b) override
    {
      assertAlways(false, "cannot compute auxiliary Jacobian vector product for AuxiliaryEquationsNone"); 
    }

  private:
    DiscPtr m_disc;
    AuxiliaryEquationsJacobiansPtr m_jacs;
};



inline AuxiliaryEquationsPtr makeAuxiliaryEquationsNone(DiscPtr disc)
{
  return std::make_shared<AuxiliaryEquationsNone>(disc);
}

#endif