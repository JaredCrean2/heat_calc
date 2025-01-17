#ifndef MESH_FACE_GROUP
#define MESH_FACE_GROUP

#include "ProjectDefs.h"
#include "mesh/tensor_product_mapper.h"

namespace Mesh {
  
// Identifies the face of an element
struct FaceSpec
{
  FaceSpec(Index el, Index el_group, LocalIndex face, LocalIndex vol_group)
    : el(el), el_group(el_group), face(face), vol_group(vol_group) {}
  Index el;             // global element index
  Index el_group;       // index within element VolumeGroup
  LocalIndex face;      // local face id
  LocalIndex vol_group; // volume group index
};


struct FaceGroup
{
  FaceGroup(const int idx, REPtr ref_el_coord, REPtr ref_el_sol,
            const TensorProductMapper& tp_mapper_coord,
            const TensorProductMapper& tp_mapper_sol,
            const std::vector<FaceSpec>& faces,
            const ArrayType<Index, 2>& nodenums,
            const std::vector<int_least8_t>& face_weights,
            bool is_dirichlet,
            bool is_boundary_surface) :
    faces(faces),
    nodenums(nodenums),
    ref_el_coord(ref_el_coord),
    ref_el_sol(ref_el_sol),
    m_idx(idx),
    m_is_dirichlet(is_dirichlet),
    m_is_boundary_surface(is_boundary_surface),
    m_tp_mapper_coord(tp_mapper_coord),
    m_tp_mapper_sol(tp_mapper_sol),
    m_face_weights(face_weights)
  {}

  std::vector<FaceSpec> faces;
  ArrayType<Index, 2> nodenums; // nfaces x npts per face
                                // TODO: now that we have the nodemaps, we
                                // don't need this large array anymore
  REPtr ref_el_coord;
  REPtr ref_el_sol;

  int getIdx() const { return m_idx; }

  bool getIsDirichlet() const { return m_is_dirichlet; }

  bool getIsBoundarySurface() const { return m_is_boundary_surface; }

  int getNumFaces() const { return faces.size();}

  int getNumCoordPtsPerFace() const { return ref_el_coord->getNumFaceNodes(); }

  int getNumSolPtsPerFace() const { return ref_el_sol->getNumFaceNodes(); }

  int_least8_t getFaceWeight(int face) const { return m_face_weights[face]; }

  // volume to face nodemap
  // nfaces_per_element x num nodes per face
  const ArrayType<LocalIndex, 2> getFaceNodesCoord() const { return ref_el_coord->getFaceNodes(); }

  const ArrayType<LocalIndex, 2> getFaceNodesSol() const { return ref_el_sol->getFaceNodes(); }

  // TensorProductMappers for the volume element this face is based on
  const TensorProductMapper& getTPMapperCoord() const { return m_tp_mapper_coord;}

  const TensorProductMapper& getTPMapperSol() const { return m_tp_mapper_sol;}

  private:
    const int m_idx;
    const bool m_is_dirichlet;
    const bool m_is_boundary_surface;
    const TensorProductMapper& m_tp_mapper_coord;
    const TensorProductMapper& m_tp_mapper_sol;
    std::vector<int_least8_t> m_face_weights;
};

}

#endif