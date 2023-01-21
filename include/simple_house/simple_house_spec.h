#ifndef SIMPLE_HOUSE_SPEC_H
#define SIMPLE_HOUSE_SPEC_H

#include <vector>
#include "mesh/mesh_generator_multi_block.h"
#include "physics/heat/HeatEquation.h"


namespace simple_house {

class SimpleHouseSpec
{
  public:
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

    Mesh::MultiBlockMeshSpec createMeshSpec();

  private:

    void addBuilding(Mesh::MultiBlockMeshSpec& mesh_spec);

    void addFoundationInsulation(Mesh::MultiBlockMeshSpec& mesh_spec);

    void addGround(Mesh::MultiBlockMeshSpec& mesh_spec);

    void setupMasks(Mesh::MultiBlockMeshSpec& mesh_spec);


};

}

#endif