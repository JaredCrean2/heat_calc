#ifndef MESH_TENSOR_PRODUCT_MAPPER_H
#define MESH_TENSOR_PRODUCT_MAPPER_H

#include "ProjectDefs.h"
#include "mesh/reference_element_interface.h"

namespace Mesh {

using REPtr = reference_element::REPtr;

// defines the mapping between regular and tensor product ordering
class TensorProductMapper
{
  public:
  // TODO: maybe ReferenceElement should not be responsible for tensor product mapping
  //       and this class should be moved into utils
    explicit TensorProductMapper(REPtr ref_el) :
      tp_nodemap(ref_el->getTPNodemap()),
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

    int getNumTPPoints() const { return m_xi.size(); }

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

}  // namespace

#endif