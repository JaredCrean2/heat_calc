#include "discretization/model.h"
#include "discretization/DirichletBC.h"

void Model::setUpBcs(DiscVectorPtr q, const Real t)
{

  for (auto& bc : m_dirichlet_bcs)
  {
    SurfDiscPtr surf = bc->getSurfDisc();
    ArrayType<Real, 1> qbc(boost::extents[surf->getNumCoordPtsPerFace()]);

    for (int face=0; face < surf->getNumFaces(); ++face)
    {
      //TODO: verify data is qbc is stored contiguously
      bc->getValue(face, t, &(qbc[0]));
      //TODO: finish this
    }


  }
}
