#include "mesh/FieldDataManager.h"
#include "discretization/disc_vector.h" 
#include "mesh/apfMDSField.h"


void FieldDataManager::attachVector(DiscVectorPtr vec, const std::string& name)
{
  // attach vector to the mesh data
  //TODO: check if vec is already attached

  // Apf won't print things unless they have the same shape as the mesh coordinate field
  apf::Field* field = apf::createPackedField(m_mesh, name.c_str(), 1, m_mesh->getShape());
  m_attached_vectors.push_back(vec);
  m_fields.push_back(field);
}


void FieldDataManager::attachNumbering(apf::ApfMDSNumbering* numbering)
{
  // attach numbering to mesh
  // Apf won't print things unless they have have the mesh coordinate field shape or
  // are element numberings
  apf::FieldShape* fshape_out = apf::getShape(numbering);
  if (apf::getShape(numbering) == m_sol_fshape || apf::getShape(numbering) == m_coord_fshape)
    fshape_out = m_mesh->getShape();

  apf::Numbering* numbering_out = apf::createNumbering(m_mesh, apf::getName(numbering), fshape_out, apf::countComponents(numbering));
  m_attached_numberings.push_back(numbering);
  m_numberings.push_back(numbering_out);
}


void FieldDataManager::syncToMesh()
{
  // sync data from attached entities to the mesh
  //prune();

  for (size_t i=0; i < m_attached_vectors.size(); ++i)
    if (auto vec = m_attached_vectors[i].lock())
    {
      copyVec(vec, m_fields[i]);
    }

  for (size_t i=0; i < m_attached_numberings.size(); ++i)
  {
    copyNumbering(m_attached_numberings[i], m_numberings[i]);
  }
}

void FieldDataManager::copyVec(DiscVectorPtr vec, apf::Field* field)
{
  if (!vec->isVectorCurrent())
    vec->syncArrayToVector();

  auto vec_vals = vec->getVector();
  apf::FieldShape* field_fshape = apf::getShape(field);

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  std::vector<double> vals(1);
  for (int dim=0; dim <= 3; ++dim)
  {
    if (!m_sol_fshape->hasNodesIn(dim))
      continue;

    it = m_mesh->begin(dim);
    while ( (e = m_mesh->iterate(it)))
    {
      if (field_fshape->countNodesOn(m_mesh->getType(e)) == 0)
        continue;

      for (int node=0; node < m_sol_fshape->countNodesOn(m_mesh->getType(e)); ++node)
      {
        if (apf::getNumber(m_is_dirichlet, e, node, 0)  )
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

void FieldDataManager::copyNumbering(apf::ApfMDSNumbering* numbering_in, apf::Numbering* numbering_out)
{
  apf::FieldShape* fshape_in = apf::getShape(numbering_in);
  apf::FieldShape* fshape_out = apf::getShape(numbering_out);

  apf::MeshEntity* e;
  apf::MeshIterator* it;
  std::vector<double> vals(apf::countComponents(numbering_out));
  for (int dim=0; dim <= 3; ++dim)
  {
    if (!fshape_in->hasNodesIn(dim) && !fshape_out->hasNodesIn(dim))
      continue;

    it = m_mesh->begin(dim);
    while ( (e = m_mesh->iterate(it)))
    {
      if (fshape_in->countNodesOn(m_mesh->getType(e)) == 0 || 
          fshape_out->countNodesOn(m_mesh->getType(e)) == 0)
        continue;

      for (int node=0; node < fshape_in->countNodesOn(m_mesh->getType(e)); ++node)
        for (int c=0; c < apf::countComponents(numbering_out); ++c)
          apf::number(numbering_out, e, node, c, apf::getNumber(numbering_in, e, node, c));
    }
    m_mesh->end(it);
  }
}