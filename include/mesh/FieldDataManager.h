#ifndef MESH_FIELD_DATA_MANAGER_H
#define MESH_FIELD_DATA_MANAGER_H

#include <apf.h>
#include <vector>
#include <memory>
#include "apfMDSField.h"

class DiscVector;
using DiscVectorPtr = std::shared_ptr<DiscVector>;

//TODO: ifdef USE_MDS_NUMBERING?

class FieldDataManager
{
  public:
    FieldDataManager(apf::Mesh* mesh, apf::FieldShape* fshape_sol, apf::ApfMDSNumbering* dof_nums, apf::ApfMDSNumbering* is_dirichlet) :
      m_mesh(mesh),
      m_sol_fshape(fshape_sol),
      m_dof_nums(dof_nums),
      m_is_dirichlet(is_dirichlet)
    {}

    void attachVector(DiscVectorPtr vec, const std::string& name);

    void attachNumbering(apf::ApfMDSNumbering* numbering);

    void syncToMesh();

  private:

    void copyVec(DiscVectorPtr vec, apf::Field* field);

    apf::Mesh* m_mesh;
    apf::FieldShape* m_sol_fshape;
    apf::ApfMDSNumbering* m_dof_nums;
    apf::ApfMDSNumbering* m_is_dirichlet;

    std::vector<std::weak_ptr<DiscVector>> m_attached_vectors;
    std::vector<apf::Field*> m_fields;
    std::vector<apf::ApfMDSNumbering*> m_attached_numberings;
    std::vector<apf::Numbering*> m_numberings;
};

#endif