#ifndef SIMPLE_HOUSE_GEOMETRY_GENERATOR_H
#define SIMPLE_HOUSE_GEOMETRY_GENERATOR_H

#include "simple_house/simple_house_spec.h"
#include "ProjectDefs.h"
#include "mesh/mesh_generator_multi_block.h"
#include "mesh/mesh_input.h"
#include "mesh/mesh.h"


namespace simple_house {

enum class SurfaceName
{
  GroundBottom,
  SouthExtWall,
  EastExtWall,
  NorthExtWall,
  WestExtWall,
  Roof,

  Floor,
  SouthIntWall,
  EastIntWall,
  NorthIntWall,
  WestIntWall,
  Ceiling,

  Lawn,

  // Non-BC surfaces
  FoundationBottom,
  FoundationInsulationBottom,
  GroundBottomBeneathFoundation
};

// The coordinate system is: North is +y, East is +x, z is into outer space
// The surfaces are 0-5 are the outer surface, 6-11 are the inner surfaces.
// The order is the same as the reference Hex elements: 0 = bottom, 5 = top,
// 1 is xz plane at y = ymin -> S
// 2 is yz plane at x = xmax -> E
// 3 is xz plane at y = ymax -> N
// 4 is yz plane at x = xmin -> W
class GeometryGenerator
{
  public:
    GeometryGenerator(SimpleHouseSpec& spec) :
      m_spec(spec.createMeshSpec())
    {
      createMeshSpec(spec);
      m_exterior_lengths = computeExteriorLengths();
      m_generator = std::make_shared<Mesh::MeshGeneratorMultiBlock>(m_spec);
      m_mesh = m_generator->generate();
      //collectGeometricFaces();
      createMeshCG();
    }

    // creates a MeshCG.  Sets up 1 volume group per mesh block.
    // Also sets up 12 face groups.  The first 6 the exterior surface, the
    // last 6 are the interior surface.  They use the standard ordering:
    // 0 = xy plane at z min
    // 1 = xz plane at y min
    // 2 = yz plane at x max
    // 3 = xz plane at y max
    // 4 = yz plane at x min
    // 5 = xy plane at z max
    std::shared_ptr<Mesh::MeshCG> getMesh() { return m_meshcg; }

    // creates the volume groups and sets the material properties
    void createVolumeGroups(std::shared_ptr<Heat::HeatEquation> heat_eqn);

    //const std::vector<int>& getExteriorGeometricFaces(int face) const { return m_exterior_geometric_faces[face]; }

    //const std::vector<int>& getInteriorGeometricFaces(int face) const { return m_interior_geometric_faces[face]; }

    Real computeInteriorVolume();

    int getSurfaceDirection(int surface);

    Real computeInteriorSurfaceArea(int direction);

    Real computeInteriorPerimeter(int direction);   

    // direction: 0 = xy plane, 1 = yz plane, 2 = xz plane
    Real computeExteriorSurfaceArea(int direction);

    Real computeExteriorPerimeter(int direction);

    Real computeLawnSurfaceArea();

    Real computeLawnPerimeter();

    int getSurfaceId(SurfaceName surf_name);

  private:
    void createMeshSpec(SimpleHouseSpec& spec);

    std::array<Real, 3> computeExteriorLengths();

    void createMeshCG();

    std::pair<Real, Real> computeLawnDimensions();

    std::shared_ptr<Mesh::MeshGeneratorMultiBlock> m_generator;
    apf::Mesh2* m_mesh;
    std::array<int, 3> m_xrange;
    std::array<int, 3> m_yrange;
    std::array<int, 3> m_zrange;
    int m_num_underground_thicknesses;
    int m_num_underground_depths;
    std::array<Real, 3> m_exterior_lengths;
    std::shared_ptr<Mesh::MeshCG> m_meshcg;

    std::vector<int> m_surface_directions;

    SimpleHouseSpec m_house_spec;
    Mesh::MultiBlockMeshSpec m_spec;
};

}  // namespace
#endif