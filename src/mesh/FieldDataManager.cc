#include "mesh/FieldDataManager.h"
#include "discretization/disc_vector.h" 


void FieldDataManager::attachVector(DiscVectorPtr vec, const std::string& name)
{
  // attach vector to the mesh data
  //TODO: check if vec is already attached

  apf::Field* field = apf::createPackedField(m_mesh, name.c_str(), 1, m_sol_fshape);
  copyVec(vec, field);
  m_attached_vectors.push_back(vec);
  m_fields.push_back(field);
}


void FieldDataManager::attachNumbering(apf::ApfMDSNumbering* numbering)
{
  // attach numbering to mesh
  apf::Numbering* numbering_out = fast_field::createSimilarNumbering(numbering);
  m_attached_numberings.push_back(numbering);
  m_numberings.push_back(numbering_out);
}


void FieldDataManager::syncToMesh()
{
  // sync data from attached entities to the mesh
  //prune();

  for (size_t i=0; i < m_attached_vectors.size(); ++i)
    if (auto vec = m_attached_vectors[i].lock())
      copyVec(vec, m_fields[i]);

  for (size_t i=0; i < m_attached_numberings.size(); ++i)
    apf::copyToApf(m_attached_numberings[i], m_numberings[i]);
}


void FieldDataManager::copyVec(DiscVectorPtr vec, apf::Field* field)
{
  if (!vec->isVectorCurrent())
    vec->syncArrayToVector();

  auto vec_vals = vec->getVector();

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  std::vector<double> vals(1);
  for (int dim=0; dim < 3; ++dim)
  {
    if (!m_sol_fshape->hasNodesIn(dim))
      continue;

    it = m_mesh->begin(dim);
    while ( (e = m_mesh->iterate(it)))
    {
      for (int node=0; node < m_sol_fshape->countNodesOn(m_mesh->getType(e)); ++node)
      {
        if (apf::getNumber(m_is_dirichlet, e, node, 0))
        {
          vals[0] = -1;
          apf::setComponents(field, e, node, vals.data());
        } else
        {
          vals[0] = vec_vals[apf::getNumber(m_dof_nums, e, node, 0)];
          apf::setComponents(field, e, node, vals.data());
        }
      }
    }
    m_mesh->end(it);
  }
}