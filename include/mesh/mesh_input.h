#ifndef MESH_INPUT_H
#define MESH_INPUT_H

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <ostream>

namespace Mesh
{


struct Point
{
  double x;
  double y;
  double z;
};

inline std::ostream& operator<<(std::ostream& os, const Point& pt)
{
  os << pt.x << ", " << pt.y << ", " << pt.z;
  return os;
}

struct ModelEntitySpec
{
  explicit ModelEntitySpec(const int dim=0, const int tag=0) :
    dim(dim),
    tag(tag)
  {}

  int dim;
  int tag;
};

inline bool operator==(const ModelEntitySpec& me1, const ModelEntitySpec& me2)
{
  return me1.dim == me2.dim && me1.tag == me2.tag;
}

inline bool operator!=(const ModelEntitySpec& me1, const ModelEntitySpec& me2)
{
  return !(me1 == me2);
}

inline std::ostream& operator<<(std::ostream& os, const ModelEntitySpec& me)
{
  os << "(" << me.dim << ", " << me.tag << ")";
  return os;
}

// Struct used for specifying the input to mesh generator
struct MeshSpec
{
  double xmin = 0;
  double xmax = 1;
  double ymin = 0;
  double ymax = 1;
  double zmin = 0;
  double zmax = 1;
  int nx = 10;
  int ny = 10;
  int nz = 10;
  int coord_order = 1;
};

inline MeshSpec getMeshSpec(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int nx, int ny, int nz, int coord_order=1)
{
  MeshSpec spec;
  spec.xmin = xmin; spec.xmax = xmax;
  spec.ymin = ymin; spec.ymax = ymax;
  spec.zmin = zmin; spec.zmax = zmax;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.coord_order = coord_order;

  return spec;
}

inline std::ostream& operator<<(std::ostream& os, const MeshSpec spec)
{
  os << "MeshSpec:" << std::endl;
  os << "  nx = " << spec.nx << ", ny = " << spec.ny << ", nz = " << spec.nz << std::endl;
  os << "  x = [" << spec.xmin << ", " << spec.xmax << "]" << std::endl;
  os << "  y = [" << spec.ymin << ", " << spec.ymax << "]" << std::endl;
  os << "  z = [" << spec.zmin << ", " << spec.zmax << "]" << std::endl;
  os << "  coord order = " << spec.coord_order;

  return os;
}



// identifies a set of geometric entities.  All entities must be
// of the same dimension
class MeshEntityGroupSpec
{
  public:
    using TStorage = std::vector<ModelEntitySpec>;

    explicit MeshEntityGroupSpec(const std::string& name) :
    m_model_entities(),
    m_is_dirichlet(false),
    m_name(name)
  {}

  void addModelEntity(const ModelEntitySpec& me, const ModelEntitySpec& parent_me = ModelEntitySpec(4, 0)) 
  {
    if (m_model_entities.size() > 0 && me.dim != m_model_entities.back().dim)
      throw std::invalid_argument("all model entities must be of same dimension");

    m_model_entities.push_back(me);
    m_parent_entities.push_back(parent_me);
  }

  void setIsDirichlet(const bool is_dirichlet)
  {
    m_is_dirichlet = is_dirichlet;
  }

  const TStorage& getModelEntities() const { return m_model_entities;}

  bool hasModelEntity(const ModelEntitySpec& me_spec) const
  {
    //TODO: store model_entities in sorted order, do binary search
    return std::find(m_model_entities.begin(), m_model_entities.end(), me_spec) != m_model_entities.end();
  }

  ModelEntitySpec getParentEntity(const ModelEntitySpec& me) const
  {
    //TODO: cache most recent index from hasModelEntity, check there first
    for (TStorage::size_type i=0; i < m_model_entities.size(); ++i)
      if (m_model_entities[i] == me)
        return m_parent_entities[i];

    throw(std::invalid_argument("requested parent model entity of model entity not on s urface"));
  }

  bool getIsDirichlet() const { return m_is_dirichlet;}

  const std::string& getName() const { return m_name;}

  int getIdx() const { return m_idx;}

  void setIdx(const int idx) { m_idx = idx;}

  private:
    TStorage m_model_entities;
    TStorage m_parent_entities;
    bool m_is_dirichlet;
    std::string m_name = "";
    int m_idx = -1;

  friend std::ostream& operator<<(std::ostream& os, const MeshEntityGroupSpec& spec);


};

inline std::ostream& operator<<(std::ostream& os, const MeshEntityGroupSpec& spec)
{
  for (int i=0; i < spec.m_model_entities.size(); ++i)
  {
    os << spec.m_model_entities[i];
    if (i <= spec.m_parent_entities.size())
      os << ", with parent " << spec.m_parent_entities[i];

    if (i < int(spec.m_model_entities.size()) - 1)
      os << std::endl;
  }
  
  return os;
}

} // namespace

#endif
