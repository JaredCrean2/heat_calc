#include "gtest/gtest.h"
#include "discretization/surface_discretization.h"
#include "physics/heat/HeatEquation.h"
#include "simple_house/geometry_generator.h"
#include "simple_house/simple_house_spec.h"
#include "mesh/mesh_input.h"

namespace {

std::array<Real, 3> computeCentroid(SurfDiscPtr surf)
{
  ArrayType<Real, 2> coords(boost::extents[surf->getNumQuadPtsPerFace()][3]);

  ArrayType<Real, 1> coords_dir(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::array<Real, 3> centroid = {0, 0, 0};
  Real area = 0;

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    surf->getFaceQuadCoords(face, coords);

    for (int dir=0; dir < 3; ++dir)
    {
      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        coords_dir[i] = coords[i][dir];

      centroid[dir] += integrateFaceScalar(surf, face, coords_dir);
    }

    for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
      coords_dir[i] = 1;

    area += integrateFaceScalar(surf, face, coords_dir);
  }

  return centroid / area;
}

std::array<Real, 3> computeAverageNormalVector(SurfDiscPtr surf)
{
  ArrayType<Real, 1> normal_dir(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::array<Real, 3> avg_normal = {0, 0, 0};

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    for (int dir=0; dir < 3; ++dir)
    {
      for (int i=0; i < surf->getNumQuadPtsPerFace(); ++i)
        normal_dir[i] = surf->normals[face][i][dir];

      avg_normal[dir] += integrateFaceScalar(surf, face, normal_dir);
    }
  }

  Real len = std::sqrt(dot(avg_normal, avg_normal));

  return avg_normal / len;
}
/*
    Mesh::MeshSpec middle_block;

    std::vector<double> horizontal_thicknesses;
    std::vector<int>    horizontal_numels;
    std::vector<Heat::VolumeGroupParams> horizontal_params;

    std::vector<double> ceiling_thicknesses;
    std::vector<int>    ceiling_numels;
    std::vector<Heat::VolumeGroupParams> ceiling_params;

    std::vector<double> foundation_thicknesses;
    std::vector<int>    foundation_numels;
    std::vector<Heat::VolumeGroupParams> foundation_params;

    std::vector<Real> foundation_insulation_thicknesses;
    std::vector<int>  foundation_insulation_numels;
    Heat::VolumeGroupParams foundation_insulation_params;

    Real ground_horizontal_thickness = 0;
    int  ground_horizontal_numel = 0;

    Real ground_depth = 0;
    int  ground_depth_numel = 0;

    Heat::VolumeGroupParams ground_params;
*/

}

TEST(SimpleHouseGeometryGenerator, Surfaces)
{

  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 2, 3, 4);
  spec.horizontal_thicknesses.push_back(5);
  spec.horizontal_numels.push_back(3);
  spec.horizontal_params.push_back({1, 1, 1});

  spec.ceiling_thicknesses.push_back(6);
  spec.ceiling_numels.push_back(3);
  spec.ceiling_params.push_back({1, 1, 1});

  spec.foundation_thicknesses.push_back(7);
  spec.foundation_numels.push_back(3);
  spec.foundation_params.push_back({1, 1, 1});

  spec.foundation_insulation_thicknesses.push_back(8);
  spec.foundation_insulation_numels.push_back(3);
  spec.foundation_insulation_params = {1, 1, 1};

  spec.ground_horizontal_thickness = 9;
  spec.ground_horizontal_numel = 3;

  spec.ground_depth = 10;
  spec.ground_depth_numel = 3;


  simple_house::GeometryGenerator generator(spec);
  std::shared_ptr<Mesh::MeshCG> mesh = generator.getMesh();
  DiscPtr disc = std::make_shared<Discretization>(mesh, 3, 3);

  std::vector<simple_house::SurfaceName> surf_names = {
                                                        simple_house::SurfaceName::GroundBottom,
                                                        simple_house::SurfaceName::SouthExtWall,
                                                        simple_house::SurfaceName::EastExtWall,
                                                        simple_house::SurfaceName::NorthExtWall,
                                                        simple_house::SurfaceName::WestExtWall,
                                                        simple_house::SurfaceName::Roof,

                                                        simple_house::SurfaceName::Floor,
                                                        simple_house::SurfaceName::SouthIntWall,
                                                        simple_house::SurfaceName::EastIntWall,
                                                        simple_house::SurfaceName::NorthIntWall,
                                                        simple_house::SurfaceName::WestIntWall,
                                                        simple_house::SurfaceName::Ceiling,

                                                        simple_house::SurfaceName::Lawn,

                                                        // Non-BC surfaces
                                                        simple_house::SurfaceName::FoundationBottom,
                                                        simple_house::SurfaceName::FoundationInsulationBottom,
                                                        simple_house::SurfaceName::GroundBottomBeneathFoundation
                                                      };

  std::vector<std::array<Real, 3>> centroids = {
                                                // exterior surfaces
                                                {1, 1.5, -25},
                                                {1, -5, 5},
                                                {7, 1.5, 5},
                                                {1, 8, 5},
                                                {-5, 1.5, 5},
                                                {1, 1.5, 10},

                                                // interior surfaces
                                                {1, 1.5, 0},
                                                {1, 0, 2},
                                                {2, 1.5, 2},
                                                {1, 3, 2},
                                                {0, 1.5, 2},
                                                {1, 1.5, 4},

                                                {1, 1.5, 0},

                                                // Non-BC surfaces
                                                {1, 1.5, -7},
                                                {1, 1.5, -15},
                                                {1, 1.5, -25}
                                              };

  std::vector<std::array<Real, 3>> normals = {
                                              { 0,  0, -1},
                                              { 0, -1,  0},
                                              { 1,  0,  0},
                                              { 0,  1,  0},
                                              {-1,  0,  0},
                                              { 0,  0,  1},

                                              { 0,  0,  1},
                                              { 0,  1,  0},
                                              {-1,  0,  0},
                                              { 0, -1,  0},
                                              { 1,  0,  0},
                                              { 0,  0, -1},

                                              { 0,  0,  1},

                                              { 0,  0, -1},
                                              { 0,  0, -1},
                                              { 0,  0, -1}
                                            };

  for (size_t i=0; i < surf_names.size(); ++i)
  {
    auto surf = disc->getSurfDisc(generator.getSurfaceId(surf_names[i]));

    auto avg_normal = computeAverageNormalVector(surf);
    auto centroid   = computeCentroid(surf);

    for (int d=0; d < 3; ++d)
      EXPECT_NEAR(avg_normal[d], normals[i][d], 1e-13);

    for (int d=0; d < 3; ++d)
      EXPECT_NEAR(centroid[d], centroids[i][d], 1e-13);      
  }
}

TEST(SimpleHouseGeometryGenerator, VolumeGroups)
{
  Real wall_kappa = 1;
  Real ceiling_kappa = 2;
  Real foundation_kappa = 3;
  Real foundation_insulation_kappa = 4;
  Real ground_kappa = 5;
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2, 0, 3, 0, 4, 2, 3, 4);
  spec.horizontal_thicknesses.push_back(5);
  spec.horizontal_numels.push_back(3);
  spec.horizontal_params.push_back({wall_kappa, 1, 1});

  spec.ceiling_thicknesses.push_back(6);
  spec.ceiling_numels.push_back(3);
  spec.ceiling_params.push_back({ceiling_kappa, 1, 1});

  spec.foundation_thicknesses.push_back(7);
  spec.foundation_numels.push_back(3);
  spec.foundation_params.push_back({foundation_kappa, 1, 1});

  spec.foundation_insulation_thicknesses.push_back(8);
  spec.foundation_insulation_numels.push_back(3);
  spec.foundation_insulation_params = {foundation_insulation_kappa, 1, 1};

  spec.ground_horizontal_thickness = 9;
  spec.ground_horizontal_numel = 3;

  spec.ground_depth = 10;
  spec.ground_depth_numel = 3;
  spec.ground_params = {ground_kappa, 1, 1};


  simple_house::GeometryGenerator generator(spec);
  std::shared_ptr<Mesh::MeshCG> mesh = generator.getMesh();
  DiscPtr disc = std::make_shared<Discretization>(mesh, 3, 3);
  std::shared_ptr<Heat::HeatEquation> heat_eqn = std::make_shared<Heat::HeatEquation>(disc);

  generator.createVolumeGroups(heat_eqn);


  auto test = [&](const std::string& name, Real kappa)
  { 
    int vol_id = mesh->getVolumeGroup(name).getIdx();
    EXPECT_EQ(heat_eqn->getVolumeGroupParams(vol_id).kappa, kappa);
  };

  // i constant
  test("volume_group-3-3-3", ground_kappa);
  test("volume_group-3-3-2", ground_kappa);
  test("volume_group-3-3-1", ground_kappa);

  test("volume_group-3-2-3", ground_kappa);
  test("volume_group-3-2-2", ground_kappa);
  test("volume_group-3-2-1", ground_kappa);

  test("volume_group-3-1-3", ground_kappa);
  test("volume_group-3-1-2", ground_kappa);
  test("volume_group-3-1-1", ground_kappa);  

  test("volume_group-30-3", ground_kappa);
  test("volume_group-30-2", ground_kappa);
  test("volume_group-30-1", ground_kappa);

  test("volume_group-31-3", ground_kappa);
  test("volume_group-31-2", ground_kappa);
  test("volume_group-31-1", ground_kappa);

  test("volume_group-32-3", ground_kappa);
  test("volume_group-32-2", ground_kappa);
  test("volume_group-32-1", ground_kappa);

  test("volume_group-33-3", ground_kappa);
  test("volume_group-33-2", ground_kappa);
  test("volume_group-33-1", ground_kappa);  

  // j constant
  test("volume_group-2-3-3", ground_kappa);
  test("volume_group-2-3-2", ground_kappa);
  test("volume_group-2-3-1", ground_kappa);

  test("volume_group-1-3-3", ground_kappa);
  test("volume_group-1-3-2", ground_kappa);
  test("volume_group-1-3-1", ground_kappa);  

  test("volume_group0-3-3", ground_kappa);
  test("volume_group0-3-2", ground_kappa);
  test("volume_group0-3-1", ground_kappa);

  test("volume_group1-3-3", ground_kappa);
  test("volume_group1-3-2", ground_kappa);
  test("volume_group1-3-1", ground_kappa);  

  test("volume_group2-3-3", ground_kappa);
  test("volume_group2-3-2", ground_kappa);
  test("volume_group2-3-1", ground_kappa);

  test("volume_group3-3-3", ground_kappa);
  test("volume_group3-3-2", ground_kappa);
  test("volume_group3-3-1", ground_kappa);

  for (int i=-3; i <= 3; ++i)
    for (int j=-3; j <= 3; ++j)   
    {
      std::string name = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "-3";
      test(name, ground_kappa);
    }


  // i == -2
  test("volume_group-2-2-3", ground_kappa);
  test("volume_group-2-2-2", foundation_insulation_kappa);
  test("volume_group-2-2-1", foundation_insulation_kappa);  

  test("volume_group-2-1-3", ground_kappa);
  test("volume_group-2-1-2", foundation_insulation_kappa);
  test("volume_group-2-1-1", foundation_insulation_kappa);

  test("volume_group-20-3", ground_kappa);
  test("volume_group-20-2", foundation_insulation_kappa);
  test("volume_group-20-1", foundation_insulation_kappa);  

  test("volume_group-21-3", ground_kappa);
  test("volume_group-21-2", foundation_insulation_kappa);
  test("volume_group-21-1", foundation_insulation_kappa);    

  test("volume_group-22-3", ground_kappa);
  test("volume_group-22-2", foundation_insulation_kappa);
  test("volume_group-22-1", foundation_insulation_kappa);    
  

  // i == 2
  test("volume_group2-2-3", ground_kappa);
  test("volume_group2-2-2", foundation_insulation_kappa);
  test("volume_group2-2-1", foundation_insulation_kappa);  

  test("volume_group2-1-3", ground_kappa);
  test("volume_group2-1-2", foundation_insulation_kappa);
  test("volume_group2-1-1", foundation_insulation_kappa);

  test("volume_group20-3", ground_kappa);
  test("volume_group20-2", foundation_insulation_kappa);
  test("volume_group20-1", foundation_insulation_kappa);  

  test("volume_group21-3", ground_kappa);
  test("volume_group21-2", foundation_insulation_kappa);
  test("volume_group21-1", foundation_insulation_kappa);    

  test("volume_group22-3", ground_kappa);
  test("volume_group22-2", foundation_insulation_kappa);
  test("volume_group22-1", foundation_insulation_kappa); 


  // j == -2
  test("volume_group-2-2-3", ground_kappa);
  test("volume_group-2-2-2", foundation_insulation_kappa);
  test("volume_group-2-2-1", foundation_insulation_kappa);

  test("volume_group-1-2-3", ground_kappa);
  test("volume_group-1-2-2", foundation_insulation_kappa);
  test("volume_group-1-2-1", foundation_insulation_kappa);  

  test("volume_group0-2-3", ground_kappa);
  test("volume_group0-2-2", foundation_insulation_kappa);
  test("volume_group0-2-1", foundation_insulation_kappa);    

  test("volume_group1-2-3", ground_kappa);
  test("volume_group1-2-2", foundation_insulation_kappa);
  test("volume_group1-2-1", foundation_insulation_kappa);  

  test("volume_group2-2-3", ground_kappa);
  test("volume_group2-2-2", foundation_insulation_kappa);
  test("volume_group2-2-1", foundation_insulation_kappa);   
  
  // j == 2
  test("volume_group-22-3", ground_kappa);
  test("volume_group-22-2", foundation_insulation_kappa);
  test("volume_group-22-1", foundation_insulation_kappa);

  test("volume_group-12-3", ground_kappa);
  test("volume_group-12-2", foundation_insulation_kappa);
  test("volume_group-12-1", foundation_insulation_kappa);  

  test("volume_group02-3", ground_kappa);
  test("volume_group02-2", foundation_insulation_kappa);
  test("volume_group02-1", foundation_insulation_kappa);    

  test("volume_group12-3", ground_kappa);
  test("volume_group12-2", foundation_insulation_kappa);
  test("volume_group12-1", foundation_insulation_kappa);  

  test("volume_group22-3", ground_kappa);
  test("volume_group22-2", foundation_insulation_kappa);
  test("volume_group22-1", foundation_insulation_kappa);      


  for (int i=-1; i <= 1; ++i)
    for (int j=-1; j <= 1; ++j)
    {
      std::string name3 = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "-3";
      std::string name2 = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "-2";
      std::string name1 = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "-1";

      test(name3, ground_kappa);
      test(name2, foundation_insulation_kappa);
      test(name1, foundation_kappa);
    }

  for (int i=-1; i <= 1; ++i)
    for (int j=-1; j <= 1; ++j)
    {
      if (i != 0 && j != 0)
      {
        std::string name = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "0";
        test(name, wall_kappa);
      }
    } 

  for (int i=-1; i <= 1; ++i)
    for (int j=-1; j <= 1; ++j)
    {
      {
        std::string name = std::string("volume_group") + std::to_string(i) + std::to_string(j) + "1";
        test(name, ceiling_kappa);
      }
    }     
}