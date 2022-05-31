#ifndef MESH_VOLUME_GROUP_H
#define MESH_VOLUME_GROUP_H

#include "ProjectDefs.h"
#include "mesh/tensor_product_mapper.h"
#include "apf.h"

namespace Mesh {

  class VolumeGroup
{
  public:
    VolumeGroup (const int idx, ArrayType<Index, 2>& nodenums, ArrayType<Real, 3>& coords,
                 const TensorProductMapper& tp_mapper_coord,
                 const TensorProductMapper& tp_mapper_sol,
                 REPtr ref_el_coord,
                 REPtr ref_el_sol,
                 std::vector<apf::MeshEntity*>& elements,
                 std::vector<int_least8_t> element_weights) :
      nodenums(nodenums),
      coords(coords),
      normals_xi(ref_el_coord->getNormals()),
      ref_el_coord(ref_el_coord),
      ref_el_sol(ref_el_sol),
      sol_degree(ref_el_sol->getDegree()),
      m_elements(elements),
      m_idx(idx),
      m_tp_mapper_coord(tp_mapper_coord),
      m_tp_mapper_sol(tp_mapper_sol),
      m_element_weights(element_weights)
    {}

    ArrayType<Index, 2> nodenums;   // nelems x npts per element
    ArrayType<Real, 3> coords;      // nelems x npts per element coord x 3
    ArrayType<Real, 2> normals_xi;  // nfaces x 3
    REPtr ref_el_coord;
    REPtr ref_el_sol;

    // const std::vector<apf::Vector3>& normals_xi;
    int sol_degree;

    std::vector<apf::MeshEntity*> m_elements;
    int getIdx() const {return m_idx;}

    int getNumElems() const { return nodenums.shape()[0];}

    int getNumSolPtsPerElement() const { return nodenums.shape()[1];}

    int getNumCoordPtsPerElement() const { return coords.shape()[1];}

    const TensorProductMapper& getTPMapperCoord() const { return m_tp_mapper_coord;}

    const TensorProductMapper& getTPMapperSol() const { return m_tp_mapper_sol;}

    int getElementWeight(int elnum) const { return m_element_weights[elnum]; }

    //TODO: is this used?
    // given derivative d/dxi, computes d/dx
    // deriv_xi: num_nodes_per_element x 3
    // deriv_x: same as above
    template <typename T>
    void rotateDerivative(ArrayType<T, 2> deriv_xi, ArrayType<T, 2> deriv_x);

  private:
    const int m_idx;
    const TensorProductMapper& m_tp_mapper_coord;
    const TensorProductMapper& m_tp_mapper_sol;
    std::vector<int_least8_t> m_element_weights;
};
}

#endif