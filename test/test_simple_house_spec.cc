#include "mesh/mesh_generator_multi_block.h"
#include "simple_house/simple_house_spec.h"
#include "gtest/gtest.h"

TEST(SimpleHouseSpec, MiddleBlockOnly)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 9.23, 0, 15.38, 0, 2.46, 5, 8, 4);

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(), 0u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 0u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(), 0u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 0u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(), 0u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 0u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(), 0u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 0u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(), 0u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 0u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(), 0u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 0u);  

  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 1);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 1);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 1);

  EXPECT_EQ(mesh_spec.create_blocks[0][0][0], false);

}


TEST(SimpleHouseSpec, OneLayerHorizontal)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2.0, 0, 3.0, 0, 4.0, 5, 8, 4);
  spec.horizontal_thicknesses.push_back(2);
  spec.horizontal_numels.push_back(5);

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 1u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 1u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(),  0u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 0u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 1u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 1u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(),  0u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 0u);

  EXPECT_EQ(mesh_spec.numel_plusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_plusy[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusy[0], 5);

  EXPECT_EQ(mesh_spec.thickness_plusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_plusy[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusy[0], 2);

  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 3);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 3);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 1);

  EXPECT_EQ(mesh_spec.create_blocks[1][1][0], false);

  for (int i=0; i <= 2; ++i)
    for (int j=0; j <= 2; ++j)
      if (i != 1 && j != 1)
        EXPECT_EQ(mesh_spec.create_blocks[i][j][0], true);
}

TEST(SimpleHouseSpec, OneLayerAllDirections)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2.0, 0, 3.0, 0, 4.0, 5, 8, 4);
  spec.horizontal_thicknesses.push_back(2);
  spec.horizontal_numels.push_back(5);

  spec.ceiling_numels.push_back(3);
  spec.ceiling_thicknesses.push_back(6);

  spec.foundation_numels.push_back(4);
  spec.foundation_thicknesses.push_back(7);

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 1u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 1u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 1u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 1u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 1u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 1u);

  EXPECT_EQ(mesh_spec.numel_plusx[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusy[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusz[0],  3);
  EXPECT_EQ(mesh_spec.numel_minusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusy[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusz[0], 4);



  EXPECT_EQ(mesh_spec.thickness_plusx[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusy[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusz[0],  6);
  EXPECT_EQ(mesh_spec.thickness_minusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusy[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusz[0], 7);


  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 3);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 3);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 3);

  EXPECT_EQ(mesh_spec.create_blocks[1][1][1], false);

  for (int i=0; i <= 2; ++i)
    for (int j=0; j <= 2; ++j)
      for (int k=0; k <= 2; ++k)
      if (i != 1 && j != 1 && k != 1)
        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], true);  
}


TEST(SimpleHouseSpec, OneLayerAllDirectionsWithFoundationInsulation)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2.0, 0, 3.0, 0, 4.0, 5, 8, 4);
  spec.horizontal_thicknesses.push_back(2);
  spec.horizontal_numels.push_back(5);

  spec.ceiling_numels.push_back(3);
  spec.ceiling_thicknesses.push_back(6);

  spec.foundation_numels.push_back(4);
  spec.foundation_thicknesses.push_back(7);

  spec.foundation_insulation_numels.push_back(6);
  spec.foundation_insulation_thicknesses.push_back(8);

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(),  2u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 2u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(),  2u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 2u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 2u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(),  2u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 2u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(),  2u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 2u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 2u);

  EXPECT_EQ(mesh_spec.numel_plusx[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusx[1],  6);

  EXPECT_EQ(mesh_spec.numel_plusy[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusy[1],  6);

  EXPECT_EQ(mesh_spec.numel_plusz[0],  3);

  EXPECT_EQ(mesh_spec.numel_minusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusx[1], 6);

  EXPECT_EQ(mesh_spec.numel_minusy[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusy[1], 6);

  EXPECT_EQ(mesh_spec.numel_minusz[0], 4);
  EXPECT_EQ(mesh_spec.numel_minusz[1], 6);


  EXPECT_EQ(mesh_spec.thickness_plusx[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusx[1],  8);

  EXPECT_EQ(mesh_spec.thickness_plusy[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusy[1],  8);

  EXPECT_EQ(mesh_spec.thickness_plusz[0],  6);

  EXPECT_EQ(mesh_spec.thickness_minusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusx[1], 8);

  EXPECT_EQ(mesh_spec.thickness_minusy[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusy[1], 8);

  EXPECT_EQ(mesh_spec.thickness_minusz[0], 7);
  EXPECT_EQ(mesh_spec.thickness_minusz[1], 8);


  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 5);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 5);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 4);

  //EXPECT_EQ(mesh_spec.create_blocks[1][1][1], false);

  for (int i=0; i <= 4; ++i)
    for (int j=0; j <= 4; ++j)
      for (int k=0; k <= 1; ++k)
      {
        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], true);
      }

  for (int i=0; i <= 4; ++i)
    for (int j=0; j <= 4; ++j)
      for (int k=2; k <= 3; ++k)
      {
        bool expected_val = true;
        if (i == 0 || i == 4 || j == 0 || j == 4)
          expected_val = false;

        if (i == 2 && j == 2 && k == 2)
          expected_val = false;

        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], expected_val);   
      }
}


TEST(SimpleHouseSpec, OneLayerAllDirectionsWithFoundationInsulationAndGroundHorizontal)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2.0, 0, 3.0, 0, 4.0, 5, 8, 4);
  spec.horizontal_thicknesses.push_back(2);
  spec.horizontal_numels.push_back(5);

  spec.ceiling_numels.push_back(3);
  spec.ceiling_thicknesses.push_back(6);

  spec.foundation_numels.push_back(4);
  spec.foundation_thicknesses.push_back(7);

  spec.foundation_insulation_numels.push_back(6);
  spec.foundation_insulation_thicknesses.push_back(8);

  spec.ground_horizontal_numel = 7;
  spec.ground_horizontal_thickness = 9;

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(),  3u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 3u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(),  3u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 3u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 2u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(),  3u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 3u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(),  3u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 3u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 2u);

  EXPECT_EQ(mesh_spec.numel_plusx[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusx[1],  6);
  EXPECT_EQ(mesh_spec.numel_plusx[2],  7);

  EXPECT_EQ(mesh_spec.numel_plusy[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusy[1],  6);
  EXPECT_EQ(mesh_spec.numel_plusy[2],  7);

  EXPECT_EQ(mesh_spec.numel_plusz[0],  3);

  EXPECT_EQ(mesh_spec.numel_minusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusx[1], 6);
  EXPECT_EQ(mesh_spec.numel_minusx[2], 7);

  EXPECT_EQ(mesh_spec.numel_minusy[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusy[1], 6);
  EXPECT_EQ(mesh_spec.numel_minusy[2], 7);

  EXPECT_EQ(mesh_spec.numel_minusz[0], 4);
  EXPECT_EQ(mesh_spec.numel_minusz[1], 6);


  EXPECT_EQ(mesh_spec.thickness_plusx[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusx[1],  8);
  EXPECT_EQ(mesh_spec.thickness_plusx[2],  9);

  EXPECT_EQ(mesh_spec.thickness_plusy[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusy[1],  8);
  EXPECT_EQ(mesh_spec.thickness_plusy[2],  9);

  EXPECT_EQ(mesh_spec.thickness_plusz[0],  6);

  EXPECT_EQ(mesh_spec.thickness_minusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusx[1], 8);
  EXPECT_EQ(mesh_spec.thickness_minusx[2], 9);

  EXPECT_EQ(mesh_spec.thickness_minusy[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusy[1], 8);
  EXPECT_EQ(mesh_spec.thickness_minusy[2], 9);

  EXPECT_EQ(mesh_spec.thickness_minusz[0], 7);
  EXPECT_EQ(mesh_spec.thickness_minusz[1], 8);


  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 7);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 7);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 4);

  //EXPECT_EQ(mesh_spec.create_blocks[1][1][1], false);

  for (int i=0; i <= 6; ++i)
    for (int j=0; j <= 6; ++j)
      for (int k=0; k <= 1; ++k)
      {
        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], true);
      }

  for (int i=0; i <= 6; ++i)
    for (int j=0; j <= 6; ++j)
      for (int k=2; k <= 3; ++k)
      {
        bool expected_val = true;
        if (i == 0 || i == 1 || i == 5 || i == 6 ||
            j == 0 || j == 1 || j == 5 || j == 6)
          expected_val = false;

        if (i == 3 && j == 3 && k == 2)
          expected_val = false;

        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], expected_val);   
      }
}

TEST(SimpleHouseSpec, OneLayerAllDirectionsWithFoundationInsulationAndGroundAll)
{
  simple_house::SimpleHouseSpec spec;
  spec.middle_block = Mesh::getMeshSpec(0, 2.0, 0, 3.0, 0, 4.0, 5, 8, 4);
  spec.horizontal_thicknesses.push_back(2);
  spec.horizontal_numels.push_back(5);

  spec.ceiling_numels.push_back(3);
  spec.ceiling_thicknesses.push_back(6);

  spec.foundation_numels.push_back(4);
  spec.foundation_thicknesses.push_back(7);

  spec.foundation_insulation_numels.push_back(6);
  spec.foundation_insulation_thicknesses.push_back(8);

  spec.ground_horizontal_numel = 7;
  spec.ground_horizontal_thickness = 9;

  spec.ground_depth_numel = 8;
  spec.ground_depth = 10;

  Mesh::MultiBlockMeshSpec mesh_spec = spec.createMeshSpec();

  EXPECT_EQ(mesh_spec.numel_plusx.size(),  3u);
  EXPECT_EQ(mesh_spec.numel_minusx.size(), 3u);
  EXPECT_EQ(mesh_spec.numel_plusy.size(),  3u);
  EXPECT_EQ(mesh_spec.numel_minusy.size(), 3u);
  EXPECT_EQ(mesh_spec.numel_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.numel_minusz.size(), 3u);

  EXPECT_EQ(mesh_spec.thickness_plusx.size(),  3u);
  EXPECT_EQ(mesh_spec.thickness_minusx.size(), 3u);
  EXPECT_EQ(mesh_spec.thickness_plusy.size(),  3u);
  EXPECT_EQ(mesh_spec.thickness_minusy.size(), 3u);
  EXPECT_EQ(mesh_spec.thickness_plusz.size(),  1u);
  EXPECT_EQ(mesh_spec.thickness_minusz.size(), 3u);

  EXPECT_EQ(mesh_spec.numel_plusx[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusx[1],  6);
  EXPECT_EQ(mesh_spec.numel_plusx[2],  7);


  EXPECT_EQ(mesh_spec.numel_plusy[0],  5);
  EXPECT_EQ(mesh_spec.numel_plusy[1],  6);
  EXPECT_EQ(mesh_spec.numel_plusy[2],  7);

  EXPECT_EQ(mesh_spec.numel_plusz[0],  3);

  EXPECT_EQ(mesh_spec.numel_minusx[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusx[1], 6);
  EXPECT_EQ(mesh_spec.numel_minusx[2], 7);

  EXPECT_EQ(mesh_spec.numel_minusy[0], 5);
  EXPECT_EQ(mesh_spec.numel_minusy[1], 6);
  EXPECT_EQ(mesh_spec.numel_minusy[2], 7);

  EXPECT_EQ(mesh_spec.numel_minusz[0], 4);
  EXPECT_EQ(mesh_spec.numel_minusz[1], 6);
  EXPECT_EQ(mesh_spec.numel_minusz[2], 8);


  EXPECT_EQ(mesh_spec.thickness_plusx[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusx[1],  8);
  EXPECT_EQ(mesh_spec.thickness_plusx[2],  9);


  EXPECT_EQ(mesh_spec.thickness_plusy[0],  2);
  EXPECT_EQ(mesh_spec.thickness_plusy[1],  8);
  EXPECT_EQ(mesh_spec.thickness_plusy[2],  9);


  EXPECT_EQ(mesh_spec.thickness_plusz[0],  6);

  EXPECT_EQ(mesh_spec.thickness_minusx[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusx[1], 8);
  EXPECT_EQ(mesh_spec.thickness_minusx[2], 9);


  EXPECT_EQ(mesh_spec.thickness_minusy[0], 2);
  EXPECT_EQ(mesh_spec.thickness_minusy[1], 8);
  EXPECT_EQ(mesh_spec.thickness_minusy[2], 9);


  EXPECT_EQ(mesh_spec.thickness_minusz[0], 7);
  EXPECT_EQ(mesh_spec.thickness_minusz[1], 8);
  EXPECT_EQ(mesh_spec.thickness_minusz[2], 10);


  EXPECT_EQ(mesh_spec.create_blocks.shape()[0], 7);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[1], 7);
  EXPECT_EQ(mesh_spec.create_blocks.shape()[2], 5);

  //EXPECT_EQ(mesh_spec.create_blocks[1][1][1], false);

  for (int i=0; i <= 6; ++i)
    for (int j=0; j <= 6; ++j)
      for (int k=0; k <= 2; ++k)
      {
        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], true);
      }

  for (int i=0; i <= 6; ++i)
    for (int j=0; j <= 6; ++j)
      for (int k=3; k <= 4; ++k)
      {
        bool expected_val = true;
        if (i == 0 || i == 1 || i == 5 || i == 6 ||
            j == 0 || j == 1 || j == 5 || j == 6)
          expected_val = false;

        if (i == 3 && j == 3 && k == 3)
          expected_val = false;

        EXPECT_EQ(mesh_spec.create_blocks[i][j][k], expected_val);   
      }
}