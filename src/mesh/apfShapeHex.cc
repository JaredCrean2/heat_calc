
#include "mesh/apfShapeHex.h"
#include "apfMesh.h"
#include <apf.h>
#include <stdexcept>
#include <string>
#include <iostream>

namespace apf {

void fail(const std::string& why)
{
  fail(why.c_str());
}

class LagrangeHexQuadratic : public FieldShape
{
  public:
    LagrangeHexQuadratic() { registerSelf(getName()); }

    const char* getName() const override {return getName_();}

    static const char* getName_() {return "LagrangeHex Quadratic";}
    
    class Quad : public EntityShape
    {
    public:
      void getValues(Mesh*, MeshEntity*,
          Vector3 const& xi, NewArray<double>& values) const
      {
        fail(std::string("getvalues not defined for ") + getName_());
      }

      void getLocalGradients(Mesh*, MeshEntity*,
          Vector3 const& xi,
          NewArray<Vector3>& grads) const
      {
        fail(std::string("getLocalGradients not defind for ") + getName_());
      }

      void getVectorValues(Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const
      {
    	  fail("getVectorValues not defined for nodal shapes");
      }

      int countNodes() const {return 9;}
    };

    class Hex : public EntityShape
    {
    public:
      void getValues(Mesh*, MeshEntity*,
          Vector3 const& xi, NewArray<double>& values) const
      {
        fail(std::string("getvalues not defined for ") + getName_());
      }

      void getLocalGradients(Mesh*, MeshEntity*,
          Vector3 const& xi,
          NewArray<Vector3>& grads) const
      {
        fail(std::string("getLocalGradients not defined for ") + getName_());
      }

      void getVectorValues(Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const
      {
    	  fail("getVectorValues not defined for nodal shapes");
      }

      int countNodes() const {return 27;}
    };

    EntityShape* getEntityShape(int type) override
    {
      static Quad quad;
      static Hex hex;
      if (type == Mesh::QUAD)
        return &quad;
      else if (type == Mesh::HEX)
        return &hex;
      //else
      throw std::runtime_error(std::string("entity type ") + apf::Mesh::typeName[type] + " not supported");
    }

    bool hasNodesIn(int dimension) override
    {
      return true;
    }
    
    int countNodesOn(int type) override
    {
      return 1;
    }

    int getOrder() override
    {
      return 2;
    }

    void getNodeXi(int type, int node, Vector3& xi) override
    {
      fail("getNodeXi not supported");
    }
};


class LagrangeHexCubic : public FieldShape
{
  public:
    LagrangeHexCubic() { registerSelf(getName()); }

    const char* getName() const override {return getName_();}

    static const char* getName_() {return "LagrangeHex Cubic";}

    
    class Quad : public EntityShape
    {
    public:
      void getValues(Mesh*, MeshEntity*,
          Vector3 const& xi, NewArray<double>& values) const
      {
        fail(std::string("getvalues not defined for ") + getName_());
      }

      void getLocalGradients(Mesh*, MeshEntity*,
          Vector3 const& xi,
          NewArray<Vector3>& grads) const
      {
        fail(std::string("getLocalGradients not defind for ") + getName_());
      }

      void getVectorValues(Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const
      {
    	  fail("getVectorValues not defined for nodal shapes");
      }

      int countNodes() const {return 16;}
    };

    class Hex : public EntityShape
    {
    public:
      void getValues(Mesh*, MeshEntity*,
          Vector3 const& xi, NewArray<double>& values) const
      {
        fail(std::string("getvalues not defined for ") + getName_());
      }

      void getLocalGradients(Mesh*, MeshEntity*,
          Vector3 const& xi,
          NewArray<Vector3>& grads) const
      {
        fail(std::string("getLocalGradients not defind for ") + getName_());
      }

      void getVectorValues(Mesh*, MeshEntity*, Vector3 const&, NewArray<Vector3>&) const
      {
    	  fail("getVectorValues not defined for nodal shapes");
      }

      int countNodes() const {return 64;}
    };

    EntityShape* getEntityShape(int type) override
    {
      static Quad quad;
      static Hex hex;
      if (type == Mesh::QUAD)
        return &quad;
      else if (type == Mesh::HEX)
        return &hex;
      else
        throw std::runtime_error(std::string("entity type ") + apf::Mesh::typeName[type] + " not supported");
    }

    bool hasNodesIn(int dimension) override
    {
      return true;
    }
    
    int countNodesOn(int type) override
    {
      if (type == Mesh::VERTEX)
        return 1;
      else if (type == Mesh::EDGE)
        return 2;
      else if (type == Mesh::QUAD)
        return 4;
      else if (type == Mesh::HEX)
        return 8;
      else
        return 0;
        //throw std::runtime_error(std::string("entity type ") + apf::Mesh::typeName[type] + " not supported");
    }

    int getOrder() override
    {
      return 3;
    }

    void getNodeXi(int type, int node, Vector3& xi) override
    {
      fail("getNodeXi not supported");
    }
};

}  // namespace

namespace Mesh {

apf::FieldShape* getLagrange(int order)
{
  static apf::LagrangeHexQuadratic quadratic;
  static apf::LagrangeHexCubic cubic;
  if (order == 1)
    return apf::getLagrange(order);
  if (order == 2)
    return &quadratic;
  if (order == 3)
    return &cubic;


  throw std::runtime_error(std::string("order ") + std::to_string(order) + " not supported");
}


}  // namespace
