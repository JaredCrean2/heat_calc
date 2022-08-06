#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include "discretization/DirichletBC.h"
#include "discretization/NeumannBC.h"
#include "discretization/source_term.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "linear_system/assembler.h"

// Evaluates a discretization for the given equation
// Produces linear systems K u = f, where K is the left-hand side
// matrix and f is the right hand side vector
class PhysicsModel
{
  public:
    explicit PhysicsModel(DiscPtr disc) :
        m_disc(disc)
    {}

    virtual ~PhysicsModel() {}

    // supports delayed initialization
    virtual void initialize() {}

    // overwrites rhs with the right hand side
    // on entry, u has the solution in vector form
    // on exit, rhs has the residual in array form
    virtual void computeRhs(DiscVectorPtr u, const Real t, DiscVectorPtr rhs) = 0;

    virtual void computeJacobian(DiscVectorPtr u, const Real t, linear_system::AssemblerPtr assembler) = 0;

    virtual void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) = 0;

    virtual void computeMassMatrix(linear_system::AssemblerPtr assembler) = 0;

    virtual DiscPtr getDiscretization() { return m_disc; }

    virtual const DiscPtr getDiscretization() const { return m_disc; }

    virtual void addDirichletBC(DirichletBCPtr bc) { m_dirichlet_bcs.push_back(bc); }

    virtual void addNeumannBC(NeumannBCPtr bc) { m_neumann_bcs.push_back(bc); }

    virtual const std::vector<DirichletBCPtr>& getDirichletBCs() const {return m_dirichlet_bcs; }

    virtual const std::vector<NeumannBCPtr>& getNeumannBCs() const {return m_neumann_bcs; }

    virtual void addSourceTerm(SourceTermPtr src) { m_source_terms.push_back(src); }

    virtual SourceTermPtr getSourceTerm(int idx)  const { return m_source_terms.at(idx); }
    
  protected:
    virtual void checkInitialization();


  private:
    DiscPtr m_disc;
    std::vector<DirichletBCPtr> m_dirichlet_bcs;
    std::vector<NeumannBCPtr> m_neumann_bcs;
    std::vector<SourceTermPtr> m_source_terms;
};

#endif