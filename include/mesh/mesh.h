#ifndef MESH_H
#define MESH_H

#include "ProjectDefs.h"
#include "mesh/mesh_input.h"
#include "mesh/reference_element.h"

#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <vector>

namespace Mesh {


bool initialize();

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


// defines the mapping between regular and tensor product ordering
class TensorProductMapper
{
  public:
  // TODO: maybe ReferenceElement should not be responsible for tensor product mapping
  //       and this class should be moved into utils
    explicit TensorProductMapper(ReferenceElement* ref_el) :
      tp_nodemap(ref_el->getTensorProductMap()),
      m_xi(ref_el->getTensorProductXi())
      //m_ref_el(ref_el)
    {}

    explicit TensorProductMapper(const std::vector<Real>& xi = {}) :
      tp_nodemap(boost::extents[xi.size()][xi.size()][xi.size()]),
      m_xi(xi)
    {
      int idx = 0;
      for (unsigned int i=0; i < xi.size(); ++i)
        for (unsigned int j=0; j < xi.size(); ++j)
          for (unsigned int k=0; k < xi.size(); ++k)
            tp_nodemap[i][j][k] = idx++;
    }

    TensorProductMapper& operator=(const TensorProductMapper& other)
    {
      tp_nodemap.resize(
        boost::extents[other.tp_nodemap.shape()[0]][other.tp_nodemap.shape()[1]][other.tp_nodemap.shape()[2]]);

      tp_nodemap = other.tp_nodemap;
      m_xi = other.m_xi;

      return *this;
    }

    const std::vector<Real> getXi() const { return m_xi;}
    
    const ArrayType<LocalIndex, 3>& getNodemap() const { return tp_nodemap;}

    boost::multi_array_types::extent_gen::gen_type<3>::type getTPShape() const
    {
      return boost::extents[tp_nodemap.shape()[0]][tp_nodemap.shape()[1]][tp_nodemap.shape()[2]];
    }

    // map 1D to 3D (tensor product) representation
    template <typename Array1D, typename Array3D>
    void mapToTP(const Array1D& from, Array3D& to) const
    {
      assert(from.num_dimensions() == 1);
      assert(to.num_dimensions()   == 3);
      for (unsigned int i=0; i < tp_nodemap.shape()[0]; ++i)
        for (unsigned int j=0; j < tp_nodemap.shape()[1]; ++j)
          for (unsigned int k=0; k < tp_nodemap.shape()[2]; ++k)
            to[i][j][k] = from[tp_nodemap[i][j][k]];
    }

    template <typename Array1D, typename Array3D>
    void mapFromTP(const Array3D& from, Array1D& to) const
    {
      assert(from.num_dimensions() == 3);
      assert(to.num_dimensions()   == 1);

      for (unsigned int i=0; i < tp_nodemap.shape()[0]; ++i)
        for (unsigned int j=0; j < tp_nodemap.shape()[1]; ++j)
          for (unsigned int k=0; k < tp_nodemap.shape()[2]; ++k)
            to[tp_nodemap[i][j][k]] = from[i][j][k];
    }

    template <typename Array3D>
    std::ostream& printTP(std::ostream& os, const Array3D& from) const
    {
      for (unsigned int i=0; i < tp_nodemap.shape()[0]; ++i)
        for (unsigned int j=0; j < tp_nodemap.shape()[1]; ++j)
        {
          for (unsigned int k=0; k < tp_nodemap.shape()[2]; ++k)
            os << from[i][j][k] << ", ";

          os << "\n";
        }

      return os;
    }



  private:
    ArrayType<LocalIndex, 3> tp_nodemap;
    std::vector<Real> m_xi;
    //ReferenceElement* m_ref_el;
};


class VolumeGroup
{
  public:
    VolumeGroup (const int idx, ArrayType<Index, 2>& nodenums, ArrayType<Real, 3>& coords,
                 const TensorProductMapper& tp_mapper_coord,
                 const TensorProductMapper& tp_mapper_sol,
                 ReferenceElement* ref_el_coord,
                 ReferenceElement* ref_el_sol,
                 std::vector<apf::MeshEntity*>& elements) :
      nodenums(nodenums),
      coords(coords),
      normals_xi(ref_el_coord->getNormals()),
      ref_el_coord(ref_el_coord),
      ref_el_sol(ref_el_sol),
      sol_degree(ref_el_sol->getDegree()),
      m_elements(elements),
      m_idx(idx),
      m_tp_mapper_coord(tp_mapper_coord),
      m_tp_mapper_sol(tp_mapper_sol)
    {}

    ArrayType<Index, 2> nodenums;   // nelems x npts per element
    ArrayType<Real, 3> coords;      // nelems x npts per element coord x 3
    ArrayType<Real, 2> normals_xi;  // nfaces x 3
    ReferenceElement* ref_el_coord;
    ReferenceElement* ref_el_sol;

    // const std::vector<apf::Vector3>& normals_xi;
    int sol_degree;

    std::vector<apf::MeshEntity*> m_elements;
    int getIdx() const {return m_idx;}

    int getNumElems() const { return nodenums.shape()[0];}

    int getNumSolPtsPerElement() const { return nodenums.shape()[1];}

    int getNumCoordPtsPerElement() const { return coords.shape()[1];}

    const TensorProductMapper& getTPMapperCoord() const { return m_tp_mapper_coord;}

    const TensorProductMapper& getTPMapperSol() const { return m_tp_mapper_sol;}

    // given derivative d/dxi, computes d/dx
    // deriv_xi: num_nodes_per_element x 3
    // deriv_x: same as above
    template <typename T>
    void rotateDerivative(ArrayType<T, 2> deriv_xi, ArrayType<T, 2> deriv_x);

  private:
    const int m_idx;
    const TensorProductMapper& m_tp_mapper_coord;
    const TensorProductMapper& m_tp_mapper_sol;
};

struct FaceGroup
{
  FaceGroup(const int idx, ReferenceElement* ref_el_coord, ReferenceElement* ref_el_sol,
            ArrayType<LocalIndex, 2> nodemap_coord,
            ArrayType<LocalIndex, 2> nodemap_sol,
            const ArrayType<LocalIndex, 2> face_tp_nodemap_coord) :
    ref_el_coord(ref_el_coord),
    ref_el_sol(ref_el_sol),
    nodemap_coord(nodemap_coord),
    nodemap_sol(nodemap_sol),
    face_tp_nodemap_coord(face_tp_nodemap_coord),
    m_idx(idx)
  {}

  std::vector<FaceSpec> faces;
  ArrayType<Index, 2> nodenums; // nfaces x npts per face
                                // TODO: now that we have the nodemaps, we
                                // don't need this large array anymore
  ReferenceElement* ref_el_coord;
  ReferenceElement* ref_el_sol;
  // volume to face nodemap
  ArrayType<LocalIndex, 2> nodemap_coord;  // nfaces_per_element x num coord pts per face
  ArrayType<LocalIndex, 2> nodemap_sol;    // nfaces_per_element x num sol pts per face

  // face tensor product nodemap
  const ArrayType<LocalIndex, 2> face_tp_nodemap_coord;

  int getIdx() const { return m_idx; }

  int getNumFaces() const { return faces.size();}

  int getNumCoordPtsPerFace() const { return nodemap_coord.shape()[1]; }

  int getNumSolPtsPerFace() const { return nodemap_sol.shape()[1];}

  private:
    const int m_idx;
};

struct ApfData
{
  explicit ApfData(apf::Mesh2* m = nullptr, apf::Numbering* dof_nums=nullptr,
                   apf::Numbering* is_dirichlet=nullptr, apf::FieldShape* sol_shape=nullptr) :
    m(m), dof_nums(dof_nums), is_dirichlet(is_dirichlet), sol_shape(sol_shape)
  {}

  apf::Mesh2* m;
  apf::Numbering* dof_nums;
  apf::Numbering* el_nums;
  apf::Numbering* is_dirichlet;  // if dof is dirichlet
  apf::FieldShape* sol_shape;    // FieldShape of solution
  apf::FieldShape* coord_shape;  // FieldShape of cooordinate field
  apf::Numbering* vol_groups;    // volume group number of elements

  //std::vector<apf::MeshEntity*> elements;
};

struct DofNumbering
{
  int sol_degree              = 0;
  int coord_degree            = 0;
  int num_dofs                = 0;
  int num_dofs_total          = 0;
  int nodes_per_element       = 0;
  int coord_nodes_per_element = 0;
  int nodes_per_face          = 0;
};

class MeshCG
{
  using SInt = std::vector<VolumeGroup>::size_type;

  public:
    MeshCG(apf::Mesh2* m,
           std::vector<MeshEntityGroupSpec> volume_group_spec,
           std::vector<MeshEntityGroupSpec> bc_spec,
           std::vector<MeshEntityGroupSpec> other_surface_spec,
           const int solution_degree, const int coord_degree);

    // getting faces
    const FaceGroup& getFaces(const MeshEntityGroupSpec& surf)
    {
      return m_all_faces.at(surf.getIdx());
    }

    const FaceGroup& getFaces(const Index idx)
    {
      return m_all_faces.at(idx);
    }

    Index getNumBCSurfaces() const {return m_bc_spec.size();}

    Index getNumSurfaces() const {return m_all_face_spec.size();}

    const MeshEntityGroupSpec& getSurface(const std::string& name) const
    {
      for (SInt idx=0; idx < m_all_face_spec.size(); ++idx)
        if (m_all_face_spec[idx].getName() == name)
          return m_all_face_spec[idx];

      throw std::invalid_argument(std::string("cannot find surface named ") + name);
    }

    // getting elements
    VolumeGroup& getElements(const MeshEntityGroupSpec& surf)
    {
      return m_vol_group.at(surf.getIdx());
    }

    VolumeGroup& getElements(const Index idx)
    {
      return m_vol_group.at(idx);
    }

    Index getNumVolumeGroups() const {return m_vol_group.size();}

    const MeshEntityGroupSpec& getVolumeGroup(const std::string& name) const
    {
      for (SInt idx=0; idx < m_volume_spec.size(); ++idx)
        if (m_all_face_spec[idx].getName() == name)
          return m_volume_spec[idx];

      throw std::invalid_argument(std::string("cannot find volume group named ") + name);
    }

    std::vector<apf::Vector3> normals_xi;

    //---------------------------------------------------------------------
    // dof functions

    // gets total number of dofs, including dirichlet BC dofs
    Index getNumTotalDofs() const {return m_dof_numbering.num_dofs_total;}

    // gets number of dofs, exluding dirichlet BC dofs
    Index getNumDofs() const {return m_dof_numbering.num_dofs;}

    // returns true if dof appears in the linear system (ie. not Dirichlet)
    bool isDofActive(const Index dof) const { return dof < getNumDofs();}

    void getElementDofs(const VolumeGroup& vol_group, int el_idx, std::vector<Index>& nodenums);

    // returns all the dofs connected to the given node (including self)
    void getDofConnectivity(const Index el, const Index node, std::vector<DofInt>);

  private:
    void setVolumeIndices();

    void setSurfaceIndices();

    void setApfData();

    void createVolumeGroups();

    void createFaceGroups();


    // input
    ApfData m_apf_data;
    DofNumbering m_dof_numbering;
    std::vector<MeshEntityGroupSpec> m_volume_spec;
    std::vector<MeshEntityGroupSpec> m_bc_spec;
    std::vector<MeshEntityGroupSpec> m_all_face_spec;

    // computed data
    std::vector<VolumeGroup> m_vol_group;
    //std::vector<std::vector<FaceGroup>> m_bc_faces;
    std::vector<FaceGroup> m_all_faces;
    std::vector<apf::MeshEntity*> m_elements;
    std::vector<Index> m_elnums_global_to_local;  // element numbers to group element numbers
    std::vector< std::vector<Index> > m_elnums_local_to_global; // (group, local element number) -> global element number
    TensorProductMapper m_tensor_product_coord_map;
    TensorProductMapper m_tensor_product_sol_map;

};


} // namespace

#endif
