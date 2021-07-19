#ifndef MODEL_H
#define MODEL_H

#include <memory>
#include "discretization/discretization.h"

class DiscVector;
class DirichletBC;
class NeumannBC;

namespace {

using DiscVectorPtr = std::shared_ptr<DiscVector>;
using DirichletBCPtr = std::shared_ptr<DirichletBC>;
using NeumannBCPtr = std::shared_ptr<NeumannBC>;

}


class Model
{
  public:

    void computeResidual(DiscVectorPtr q, DiscVectorPtr res);

    void computeJacobian();

  private:

    void setUpBcs(DiscVectorPtr q, const Real t);

    DiscPtr m_disc;
    std::vector<DirichletBCPtr> m_dirichlet_bcs;
    std::vector<NeumannBCPtr> m_neumann_bcs;
};

#endif



