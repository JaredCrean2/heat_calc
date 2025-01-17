#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include "discretization/DirichletBC.h"
#include "discretization/NeumannBC.h"
#include "discretization/source_term.h"
#include "discretization/disc_vector.h"
#include "discretization/discretization.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix.h"
#include "AuxiliaryEquations.h"
#include "post_processor_manager.h"

enum class JacobianTerms
{
  Linear,
  Nonlinear,
  All
};

// Evaluates a discretization for the given equation
// Produces linear systems K u = f, where K is the left-hand side
// matrix and f is the right hand side vector
class PhysicsModel
{
  public:
    explicit PhysicsModel(DiscPtr disc, MPI_Comm comm=MPI_COMM_WORLD) :
        m_disc(disc),
        m_aux_equations_none(makeAuxiliaryEquationsNone(disc)),
        m_comm(comm)

    {}

    virtual ~PhysicsModel() {}

    // supports delayed initialization
    virtual void initialize() {}

    // overwrites rhs with the right hand side
    // on entry, u has the solution in vector form
    // on exit, rhs has the residual in array form
    virtual void computeRhs(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, DiscVectorPtr rhs) = 0;

    virtual void computeJacobian(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, const Real t, linear_system::AssemblerPtr assembler,
                                 JacobianTerms terms=JacobianTerms::All) = 0;

    virtual void applyMassMatrix(DiscVectorPtr vec_in, DiscVectorPtr vec_out) = 0;

    virtual void computeMassMatrix(linear_system::AssemblerPtr assembler) = 0;

    virtual DiscPtr getDiscretization() { return m_disc; }

    virtual const DiscPtr getDiscretization() const { return m_disc; }

    virtual void addDirichletBC(DirichletBCPtr bc) { m_dirichlet_bcs.push_back(bc); }

    virtual void addNeumannBC(NeumannBCPtr bc) { m_neumann_bcs.push_back(bc); }

    virtual const std::vector<DirichletBCPtr>& getDirichletBCs() const {return m_dirichlet_bcs; }

    virtual const std::vector<NeumannBCPtr>& getNeumannBCs() const {return m_neumann_bcs; }

    virtual void addSourceTerm(int vol_disc, SourceTermPtr src) { m_source_terms[vol_disc] = src; }

    virtual SourceTermPtr getSourceTerm(int idx)  const { return m_source_terms.at(idx); }

    bool hasSourceTerm(size_t vol_disc) const { return m_source_terms.count(vol_disc) > 0; }

    virtual AuxiliaryEquationsPtr getAuxEquations() { return m_aux_equations_none; }

    void setPostProcessors(physics::PostProcessorManagerPtr postprocs) { m_postprocessors = postprocs; }

    physics::PostProcessorManagerPtr getPostProcessors() { return m_postprocessors; }

    void runPostProcessors(int timestep, DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t)
    {
      if (m_postprocessors)
        m_postprocessors->runPostProcessors(timestep, u, u_aux, t);
    }

    MPI_Comm getComm() const { return m_comm; }

  protected:
    virtual void checkInitialization();


  private:
    DiscPtr m_disc;
    std::vector<DirichletBCPtr> m_dirichlet_bcs;
    std::vector<NeumannBCPtr> m_neumann_bcs;
    std::map<int, SourceTermPtr> m_source_terms;
    AuxiliaryEquationsPtr m_aux_equations_none;
    physics::PostProcessorManagerPtr m_postprocessors;
    MPI_Comm m_comm;
};

#endif