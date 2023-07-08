#ifndef REFERENCE_ELEMENT_APF_H
#define REFERENCE_ELEMENT_APF_H

#include "error_handling.h"
#include "mesh/reference_element_interface.h"

#include "apfMesh.h"
#include "apf.h"
#include "apfShape.h"
#include <stdexcept>
#include <string>
#include <iostream>

namespace apf {

//TODO: get rid of other definition of fail, and then remove the Impl namespace
namespace Impl {
void fail(const std::string& why);
}

// make an apf::EntityShape from a ReferenceElement
class EntityShapeRefEl : public EntityShape
{
  public:
    EntityShapeRefEl(std::shared_ptr<reference_element::ReferenceElement> ref_el, int dim, int entity) :
      m_ref_el(ref_el),
      m_dim(dim)
    {
      auto ge = ref_el->getEntity(dim, entity);
      
      m_nnodes = ref_el->getNumNodes(m_dim);
      for (int d=0; d < m_dim; ++d)
      {
        auto entities_i = reference_element::getDownward(ge, d);
        m_nnodes += entities_i.size() * ref_el->getNumNodes(d);
      }
    }

    void getValues(Mesh*, MeshEntity*,
        Vector3 const& xi, NewArray<double>& values) const
    {
      Impl::fail(std::string("getValues not defined"));
    }

    void getLocalGradients(Mesh*, MeshEntity*,
        Vector3 const& xi,
        NewArray<Vector3>& grads) const
    {
      Impl::fail(std::string("getLocalGradients not defined"));
    }

    void getVectorValues(Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const
    {
      fail("getVectorValues not defined for nodal shapes");
    }

    int countNodes() const {return m_nnodes;}

    void alignSharedNodes(apf::Mesh* m, MeshEntity*elem, MeshEntity* boundary, int order[])
    {
      assertAlways(apf::Mesh::typeDimension[m->getType(boundary)] == m_dim, "the boundary entity should be the same dimension as the EntityShape");
      int which, rotate;
      bool flip;
      apf::getAlignment(m, elem, boundary, which, flip, rotate);
      
      int nodes_per_direction = m_ref_el->getNumNodes(m_dim);
      if (m_dim == 0)
      {
        assertAlways(!flip && rotate == 0, "flip and rotate must be zero for vertices");
        assertAlways(nodes_per_direction >= 0 && nodes_per_direction <= 1, "There can be either zero or one nodes per direction on vertices");
        if (nodes_per_direction == 1)
          order[0] = 0;
      } else if (m_dim == 1)
      {
        assertAlways(rotate == 0, "rotate must be zero for edges");
        if (flip)
          for (int i=0; i < nodes_per_direction; ++i)
            order[i] = nodes_per_direction - i - 1;
        else
          for (int i=0; i < nodes_per_direction; ++i)
            order[i] = i;
      } else  // face
      {
        for (int node=0; node < nodes_per_direction; ++node)
        {
          // Note: this is based on how the face nodes are defined in tensor product order
          const int i = node % nodes_per_direction;
          const int j = node / nodes_per_direction;
          int i2 = -1, j2 = -1;
          if (!flip)
          {
            switch (rotate)
            {
              case 0: {i2 = i; j2 = j; break;}
              case 1: {i2 = j; j2 = nodes_per_direction - i - 1; break;}
              case 2: {i2 = nodes_per_direction - i - 1; j2 = nodes_per_direction - j - 1; break;}
              case 3: {i2 = nodes_per_direction - j - 1; j2 = i; break; }
              default:
                throw std::runtime_error("rotate value not recognized");
            };
          } else
          {
            switch (rotate)
            {
              case 0: {i2 = i; j2 = nodes_per_direction -j - 1; break;}
              case 1: {i2 = nodes_per_direction - j - 1; j2 = nodes_per_direction - i - 1; break;}
              case 2: {i2 = nodes_per_direction - i - 1; j2 = j; break;}
              case 3: {i2 = j; j2 = i; break; }
              default:
                throw std::runtime_error("rotate value not recognized");
            };
          }

          assert(i2 >= 0 && j2 >= 0);
          order[node] = i2 + nodes_per_direction * j2;
        }
      }
    }

  private:
    std::shared_ptr<reference_element::ReferenceElement> m_ref_el;
    int m_dim;
    int m_nnodes;
};



class FieldShapeRefEl: public FieldShape
{
  public:
    FieldShapeRefEl(const std::string& name, std::shared_ptr<reference_element::ReferenceElement> ref_el,
                         std::array<std::shared_ptr<EntityShape>, 4> entityshapes) :
    m_name(name),
    m_ref_el(ref_el),
    m_entityshapes(entityshapes)
    { 
      registerSelf(getName());
    }

    const char* getName() const override {return m_name.c_str(); }    

    EntityShape* getEntityShape(int type) override
    {
      if (type == Mesh::VERTEX)
        return m_entityshapes[0].get();
      else if (type == Mesh::EDGE)
        return m_entityshapes[1].get();
      else if (type == Mesh::QUAD)
        return m_entityshapes[2].get();
      else if (type == Mesh::HEX)
        return m_entityshapes[3].get();
      //else
      throw std::runtime_error(std::string("entity type ") + apf::Mesh::typeName[type] + " not supported");
    }

    bool hasNodesIn(int dimension) override
    {
      return m_ref_el->getNumNodes(dimension) > 0;
    }
    
    int countNodesOn(int type) override
    {
      // apf sometimes queries this function for all types, even of none of those exist
      // in the mesh
      if (type == Mesh::VERTEX || type == Mesh::EDGE || type == Mesh::QUAD || type == Mesh::HEX)
      {
        int dim = apf::Mesh::typeDimension[type];
        return m_ref_el->getNumNodes(dim);
      } else
      {
        return 0;
      }
    }

    int getOrder() override
    {
      return m_ref_el->getNumNodes(1) + 1;
    }

    void getNodeXi(int type, int node, Vector3& xi) override
    {
      Impl::fail("getNodeXi not supported");
    }

  private:
    std::string m_name;
    std::shared_ptr<reference_element::ReferenceElement> m_ref_el;
    std::array<std::shared_ptr<EntityShape>, 4> m_entityshapes;
};


std::shared_ptr<FieldShapeRefEl> getHexFieldShape(std::shared_ptr<reference_element::ReferenceElement> ref_el);

}  // namespace
#endif