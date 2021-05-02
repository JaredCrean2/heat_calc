#ifndef MESH_INPUT_H
#define MESH_INPUT_H

#include <vector>

namespace Mesh
{

struct Point
{
  double x;
  double y;
  double z;
};

struct ModelEntitySpec
{
  explicit ModelEntitySpec(const int dim=0, const int tag=0) :
    dim(dim),
    tag(tag)
  {}

  int dim;
  int tag;
};

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



// identifies a set of geometric entities.  All entities must be
// of the same dimension
class Surface
{
  public:
    using TStorage = std::vector<ModelEntitySpec>;
  explicit Surface(const std::string& name) :
    m_model_entities(),
    m_is_dirichlet(false),
    m_name(name)
  {}

  void addModelEntity(const ModelEntitySpec& me) 
  {
    if (m_model_entities.size() > 0 && me.dim != m_model_entities.back().dim)
      throw std::invalid_argument("all model entities must be of same dimension");

    m_model_entities.push_back(me);
  }

  void setIsDirichlet(const bool is_dirichlet)
  {
    m_is_dirichlet = is_dirichlet;
  }

  const TStorage& getModelEntities() const { return m_model_entities;}

  bool getIsDirichlet() const { return m_is_dirichlet;}

  const std::string& getName() const { return m_name;}

  int getIdx() const { return m_idx;}

  private:
    std::vector<ModelEntitySpec> m_model_entities;
    bool m_is_dirichlet;
    std::string m_name = "";
    int m_idx = -1;
};

} // namespace

#endif
