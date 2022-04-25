#ifndef APF_MDS_FIELD_H
#define APF_MDS_FIELD_H

#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfMDS.h"
#include "mds_apf.h"
#include "apfShape.h"
#include "apfNumbering.h"
#include <stdexcept>
#include <vector>
#include <cassert>
#include <iostream>

// a faster implementation of apf::Fields/apf::Numbering specific to the
// MDS data structure

namespace fast_field {

//TODO: call reorderMdsMesh to make sure the MDS indexing is correct

class MDSMeshWrapper
{
  public:
    explicit MDSMeshWrapper(apf::Mesh2* m):
      mesh(m)
    {}

    using LargeInt = mds_id; //TODO: mds_id?

    apf::Mesh2* mesh;

    // gets an upper bounds on the value returned by getEntityIndex()
    std::vector<LargeInt> getEntityIndexUpperBounds()
    {
      const mds_apf* mds_apf = apf::getMDS(mesh);
      std::vector<LargeInt> bounds;
      for (int t=0; t < MDS_TYPES; ++t)
        bounds.push_back(mds_apf->mds.cap[t]);

      return bounds;
    }

    // return upper bounds on the value returned by getEntityType
    int getEntityTypeUpperBound()
    {
      return MDS_TYPES;
    }

    int mapType(apf::Mesh::Type type)
    {
      return apf2mds(type);
    }

    // get the type of the entity
    int getEntityType(apf::MeshEntity* e)
    {
      mds_id id = apf::fromEnt(e);
      return mds_type(id);
    }

    // get the type-specific index of the entity
    LargeInt getEntityIndex(apf::MeshEntity* e)
    {
      mds_id id = apf::fromEnt(e);
      return mds_index(id);
    }
};


template <typename T, typename MeshWrapper>
class ApfMDSField
{
  public:
    ApfMDSField(apf::Mesh2* m, const std::string& name, apf::FieldShape* fshape, int ncomp) :
      m_meshwrapper(m),
      m_name(name),
      m_fshape(fshape),
      m_ncomp(ncomp)
    {
      int ntypes = m_meshwrapper.getEntityTypeUpperBound();
      auto counts = m_meshwrapper.getEntityIndexUpperBounds();
      m_nodes_per_type.resize(ntypes);
      m_data.resize(ntypes);
      for (int i=0; i < apf::Mesh::TYPES; ++i)
      {
        apf::Mesh::Type apf_type = static_cast<apf::Mesh::Type>(i);
        int mds_type = m_meshwrapper.mapType(apf_type);
        int nodes = fshape->countNodesOn(i);
        m_nodes_per_type[mds_type] = nodes;
        m_data[mds_type].resize(m_ncomp * nodes * counts[mds_type]);
      }
    }

    virtual ~ApfMDSField() {}

    apf::Mesh* getMesh() {return m_meshwrapper.mesh; }

    const std::string& getName() const { return m_name;}

    apf::FieldShape* getFieldShape() { return m_fshape; }

    int getNumComponents() const { return m_ncomp; }

    T& operator()(apf::MeshEntity* e, int node, int component)
    {
      int mds_type = m_meshwrapper.getEntityType(e);
      int nodes_per_type = m_nodes_per_type[mds_type];
      assert(node >= 0 && node < nodes_per_type);
      assert(component >=0 && component < m_ncomp);

      LargeInt idx  = m_meshwrapper.getEntityIndex(e) * nodes_per_type * m_ncomp + node * m_ncomp + component;
      assert(idx >= 0 && idx < m_data[mds_type].size());
#ifdef NDEBUG
      return m_data[mds_type][idx];
#else
      return m_data[mds_type].at(idx);
#endif
    }

    void set(const T& val)
    {
      for (auto& arr : m_data)
        for (auto& c : arr)
          c = val;
    }

  protected:
    using LargeInt = typename MeshWrapper::LargeInt;

    MeshWrapper m_meshwrapper;
    std::string m_name;
    apf::FieldShape* m_fshape;
    int m_ncomp;
    // Note: these arrays are indexed by *mds* types
    std::vector<int> m_nodes_per_type;
    std::vector<std::vector<T>> m_data;
};

template <typename MeshWrapper>
class ApfMDSNumberingSpec : public ApfMDSField<int, MeshWrapper>
{
  using Base = ApfMDSField<int, MeshWrapper>;
  using Base::set;
  using Base::m_data;
  public:
    ApfMDSNumberingSpec(apf::Mesh2* m, const std::string& name, apf::FieldShape* fshape, int ncomp) :
      Base(m, name, fshape, ncomp)
    {
      setUnnumbered();
    }

    void fix(apf::MeshEntity* e, int node, int component, bool fixed)
    {
      if (fixed)
        (*this)(e, node, component) = m_fixed;
    }

    bool isFixed(apf::MeshEntity* e, int node, int component)
    {
      return (*this)(e, node, component) == m_fixed;
    }

    bool isNumbered(apf::MeshEntity* e, int node, int component)
    {
      return (*this)(e, node, component) > 0;  // not fixed or unnumbered
    }

    int countFixed()
    {
      int count;
      for (auto& arr : m_data)
        for (auto& c : arr)
          if (c == m_fixed)
            ++count;
      return count;
    }

    void setUnnumbered()
    {
      set(m_unnumbered);
    }

  private:
    int m_unnumbered = -1;
    int m_fixed      = -2;
};


}  // namespace

namespace apf {



using ApfMDSNumbering = fast_field::ApfMDSNumberingSpec<fast_field::MDSMeshWrapper>;

inline ApfMDSNumbering* createNumberingMDS(Mesh2* mesh, const char* name, FieldShape* shape, int component)
{
  return new ApfMDSNumbering(mesh, name, shape, component);
}

inline void destroyNumbering(ApfMDSNumbering* n)
{
  delete n;
}

inline void fix(ApfMDSNumbering* n, MeshEntity*e, int node, int component, bool fixed)
{
  n->fix(e, node, component, fixed);
}

inline bool isFixed(ApfMDSNumbering* n, MeshEntity*e, int node, int component)
{
  return n->isFixed(e, node, component);
}

inline bool isNumbered(ApfMDSNumbering* n, MeshEntity*e, int node, int component)
{
  return n->isNumbered(e, node, component);
}

inline void number(ApfMDSNumbering* n, MeshEntity*e, int node, int component, int number)
{
  (*n)(e, node, component) = number;
}

inline int getNumber(ApfMDSNumbering* n, MeshEntity*e, int node, int component)
{
  return (*n)(e, node, component);
}

inline FieldShape* getShape(ApfMDSNumbering* n)
{
  return n->getFieldShape();
}

inline const char* getName(ApfMDSNumbering* n)
{
  return n->getName().c_str();
}

inline Mesh* getMesh(ApfMDSNumbering* n)
{
  return n->getMesh();
}

inline int countComponents(ApfMDSNumbering* n)
{
  return n->getNumComponents();
}

}

#endif