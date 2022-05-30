#include "gtest/gtest.h"
#include "test_helper.h"
#include "mesh_helper.h"
#include "discretization/disc_vector.h"


namespace {
  class DiscVectorTester : public ::testing::Test,
                           public StandardDiscSetup
{
  protected:
    DiscVectorTester()
    {
      setup();
      disc_vec = makeDiscVector(disc);
    }

    DiscVectorPtr disc_vec;
};

Real func(const Real x, const Real y, const Real z)
{
  return x + 2*y + 3*z;
}

}  // namespace


TEST_F(DiscVectorTester, Sizes)
{
  EXPECT_EQ(disc_vec->getNumArrays(), 1u);
}

TEST_F(DiscVectorTester, Accesors)
{
  EXPECT_NO_THROW(disc_vec->getVector());
  EXPECT_NO_THROW(disc_vec->getArray(0));

  disc_vec->markArrayModified();
  EXPECT_NO_THROW(disc_vec->getArray(0));
  EXPECT_ANY_THROW(disc_vec->getVector());

  disc_vec = makeDiscVector(disc);
  disc_vec->markVectorModified();
  EXPECT_NO_THROW(disc_vec->getVector());
  EXPECT_ANY_THROW(disc_vec->getArray(0));


  disc_vec = makeDiscVector(disc);
  for (int i=1; i < 5; ++i)
    EXPECT_ANY_THROW(disc_vec->getArray(i));

  for (int i=-5; i < 0; ++i)
    EXPECT_ANY_THROW(disc_vec->getArray(i));
}


TEST_F(DiscVectorTester, AccesorsConst)
{
  {
    const auto disc_vec2 = disc_vec;
    EXPECT_NO_THROW(disc_vec2->getVector());
    EXPECT_NO_THROW(disc_vec2->getArray(0));

    disc_vec->markArrayModified();
    EXPECT_NO_THROW(disc_vec2->getArray(0));
    EXPECT_ANY_THROW(disc_vec2->getVector());
  }

  {
    const auto disc_vec2 = makeDiscVector(disc);
    disc_vec2->markVectorModified();
    EXPECT_NO_THROW(disc_vec2->getVector());
    EXPECT_ANY_THROW(disc_vec2->getArray(0));
  }


  {
    const auto disc_vec2 = makeDiscVector(disc);
    for (int i=1; i < 5; ++i)
      EXPECT_ANY_THROW(disc_vec2->getArray(i));

    for (int i=-5; i < 0; ++i)
      EXPECT_ANY_THROW(disc_vec2->getArray(i));
  }
}

TEST_F(DiscVectorTester, SetConstant)
{
  disc_vec->set(1);
  auto& vec = disc_vec->getVector();
  for (unsigned int i=0; i < vec.shape()[0]; ++i)
    EXPECT_EQ(vec[i], 1.0);

  auto& arr = disc_vec->getArray(0);
  for (unsigned int i=0; i < arr.shape()[0]; ++i)
    for (unsigned int j=0; j < arr.shape()[1]; ++j)
      EXPECT_EQ(arr[i][j], 1.0);
}

TEST_F(DiscVectorTester, SetFunc)
{
  disc_vec->set(1.0);
  disc_vec->setFunc(func);
  auto dof_numbering = disc->getDofNumbering();

  auto& dofs = dof_numbering->getDofs(0);
  auto& coords = disc->getVolDisc(0)->vol_group.coords;

  auto& vec = disc_vec->getVector();
  auto& arr = disc_vec->getArray(0);
  for (unsigned int i=0; i < dofs.shape()[0]; ++i)
    for (unsigned int j=0; j < dofs.shape()[1]; ++j)
    {
      Index dof = dofs[i][j];
      Real val_ex = func(coords[i][j][0], coords[i][j][1], coords[i][j][2]);

      if (dof_numbering->isDofActive(dof))
        EXPECT_EQ(vec[dof], val_ex);

      EXPECT_EQ(arr[i][j], val_ex);
    }
}

TEST_F(DiscVectorTester, syncArrayToVector)
{
  disc_vec->setFunc(func);
  auto& vec = disc_vec->getVector();
  for (unsigned int i=0; i < vec.shape()[0]; ++i)
    vec[i] = 0;

  disc_vec->markArrayModified();
  disc_vec->syncArrayToVector(Assign2<Real>());

  auto dof_numbering = disc->getDofNumbering();
  auto& dofs = dof_numbering->getDofs(0);
  auto& coords = disc->getVolDisc(0)->vol_group.coords;

  for (unsigned int i=0; i < dofs.shape()[0]; ++i)
    for (unsigned int j=0; j < dofs.shape()[1]; ++j)
    {
      Index dof = dofs[i][j];
      Real val_ex = func(coords[i][j][0], coords[i][j][1], coords[i][j][2]);

      if (dof_numbering->isDofActive(dof))
        EXPECT_EQ(vec[dof], val_ex);
    }
}

TEST_F(DiscVectorTester, syncVectorToArray)
{  
  auto dof_numbering = disc->getDofNumbering();
  auto& dofs = dof_numbering->getDofs(0);
  auto& coords = disc->getVolDisc(0)->vol_group.coords;

  disc_vec->setFunc(func);
  auto& arr = disc_vec->getArray(0);
  
  for (unsigned int i=0; i < dofs.shape()[0]; ++i)
    for (unsigned int j=0; j < dofs.shape()[1]; ++j)
      arr[i][j] = 0;

  disc_vec->markVectorModified();
  disc_vec->syncVectorToArray();

  for (unsigned int i=0; i < dofs.shape()[0]; ++i)
    for (unsigned int j=0; j < dofs.shape()[1]; ++j)
    {
      Index dof = dofs[i][j];
      Real val_ex = func(coords[i][j][0], coords[i][j][1], coords[i][j][2]);
      if (!(dof_numbering->isDofActive(dof)))
        val_ex = 0;

      EXPECT_EQ(arr[i][j], val_ex);        
    }
}