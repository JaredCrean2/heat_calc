#include "gtest/gtest.h"
#include "physics/heat/environment_interface.h"

TEST(EnvironmentInterface, Addition)
{
  Heat::EnvironmentData a, b, c;

  a.air_temp = 2;
  a.air_speed = 3;
  a.air_direction = {1, 2, 3};
  a.ir_horizontal_radiation = 4;
  a.direct_normal_radiation = 5;
  a.diffuse_radiation = 6;

  b.air_temp = 7;
  b.air_speed = 8;
  b.air_direction = {4, 5, 6};
  b.ir_horizontal_radiation = 9;
  b.direct_normal_radiation = 10;
  b.diffuse_radiation = 11;

  c = a + b;

  EXPECT_EQ(c.air_temp, 9);
  EXPECT_EQ(c.air_speed, 11);
  EXPECT_EQ(c.air_direction[0], 5);
  EXPECT_EQ(c.air_direction[1], 7);
  EXPECT_EQ(c.air_direction[2], 9);
  EXPECT_EQ(c.ir_horizontal_radiation, 13);
  EXPECT_EQ(c.direct_normal_radiation, 15);
  EXPECT_EQ(c.diffuse_radiation, 17);
}

TEST(EnvironmentInterface, Subtraction)
{
  Heat::EnvironmentData a, b, c;

  a.air_temp = 2;
  a.air_speed = 3;
  a.air_direction = {1, 2, 3};
  a.ir_horizontal_radiation = 4;
  a.direct_normal_radiation = 5;
  a.diffuse_radiation = 6;

  b.air_temp = 7;
  b.air_speed = 8;
  b.air_direction = {4, 5, 6};
  b.ir_horizontal_radiation = 9;
  b.direct_normal_radiation = 10;
  b.diffuse_radiation = 11;

  c = a - b;

  EXPECT_EQ(c.air_temp, -5);
  EXPECT_EQ(c.air_speed, -5);
  EXPECT_EQ(c.air_direction[0], -3);
  EXPECT_EQ(c.air_direction[1], -3);
  EXPECT_EQ(c.air_direction[2], -3);
  EXPECT_EQ(c.ir_horizontal_radiation, -5);
  EXPECT_EQ(c.direct_normal_radiation, -5);
  EXPECT_EQ(c.diffuse_radiation, -5);
}

TEST(EnvironmentInterface, Multiplication)
{
  Heat::EnvironmentData a, c, d;
  Real b = 2;

  a.air_temp = 2;
  a.air_speed = 3;
  a.air_direction = {1, 2, 3};
  a.ir_horizontal_radiation = 4;
  a.direct_normal_radiation = 5;
  a.diffuse_radiation = 6;

  c = a * b;

  EXPECT_EQ(c.air_temp, a.air_temp*b);
  EXPECT_EQ(c.air_speed, a.air_speed*b);
  EXPECT_EQ(c.air_direction[0], a.air_direction[0]*b);
  EXPECT_EQ(c.air_direction[1], a.air_direction[1]*b);
  EXPECT_EQ(c.air_direction[2], a.air_direction[2]*b);
  EXPECT_EQ(c.ir_horizontal_radiation, a.ir_horizontal_radiation*b);
  EXPECT_EQ(c.direct_normal_radiation, a.direct_normal_radiation*b);
  EXPECT_EQ(c.diffuse_radiation, a.diffuse_radiation*b);

  d = b * a;

  EXPECT_EQ(d.air_temp, a.air_temp*b);
  EXPECT_EQ(d.air_speed, a.air_speed*b);
  EXPECT_EQ(d.air_direction[0], a.air_direction[0]*b);
  EXPECT_EQ(d.air_direction[1], a.air_direction[1]*b);
  EXPECT_EQ(d.air_direction[2], a.air_direction[2]*b);
  EXPECT_EQ(d.ir_horizontal_radiation, a.ir_horizontal_radiation*b);
  EXPECT_EQ(d.direct_normal_radiation, a.direct_normal_radiation*b);
  EXPECT_EQ(d.diffuse_radiation, a.diffuse_radiation*b);  
}

TEST(EnvironmentInterface, Division)
{
  Heat::EnvironmentData a, c;
  Real b = 2;

  a.air_temp = 2;
  a.air_speed = 3;
  a.air_direction = {1, 2, 3};
  a.ir_horizontal_radiation = 4;
  a.direct_normal_radiation = 5;
  a.diffuse_radiation = 6;

  c = a / (1/b);

  EXPECT_EQ(c.air_temp, a.air_temp*b);
  EXPECT_EQ(c.air_speed, a.air_speed*b);
  EXPECT_EQ(c.air_direction[0], a.air_direction[0]*b);
  EXPECT_EQ(c.air_direction[1], a.air_direction[1]*b);
  EXPECT_EQ(c.air_direction[2], a.air_direction[2]*b);
  EXPECT_EQ(c.ir_horizontal_radiation, a.ir_horizontal_radiation*b);
  EXPECT_EQ(c.direct_normal_radiation, a.direct_normal_radiation*b);
  EXPECT_EQ(c.diffuse_radiation, a.diffuse_radiation*b);
}