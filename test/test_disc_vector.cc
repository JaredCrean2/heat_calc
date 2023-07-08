#include "gtest/gtest.h"
#include <apfMesh.h>
#include "discretization/DirichletBC_defs.h"
#include "test_helper.h"
#include "mesh_helper.h"
#include "discretization/disc_vector.h"
#include "discretization/DirichletBC.h"


namespace {

Real func(const Real x, const Real y, const Real z)
{
  return x + 2*y + 3*z;
}

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

class DiscVectorDirichletTester : public ::testing::Test,
                                  public StandardDiscSetupMulti
{
  protected:
    DiscVectorDirichletTester()
    {
      std::vector<bool> is_dirichlet(11, false);
      is_dirichlet[1] = true;

      setup(3, 1, getStandardMeshSpecs(), is_dirichlet); // setup(3, 3, getStandardMeshSpec(), is_dirichlet);
      disc_vec = makeDiscVector(disc);

      auto f = [](Real x, Real y, Real z, Real t) { return func(x, y, z); };
      m_bc = makeDirichletBCMMS(surf_discs[1], f);
    }

    DiscVectorPtr disc_vec;
    DirichletBCPtr m_bc;
};


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

TEST_F(DiscVectorDirichletTester, syncVectorToArrayDirichlet)
{  
  auto dof_numbering = disc->getDofNumbering();
  //auto& coords = disc->getVolDisc(0)->vol_group.coords;

  disc_vec->setFunc(func);
  for (int vol_block=0; vol_block < 2; ++vol_block)
  {
    auto& dofs = dof_numbering->getDofs(vol_block);
    auto& arr = disc_vec->getArray(vol_block);
    
    for (unsigned int i=0; i < dofs.shape()[0]; ++i)
      for (unsigned int j=0; j < dofs.shape()[1]; ++j)
        arr[i][j] = 0;
  }

  disc_vec->markVectorModified();
  disc_vec->syncVectorToArray();
  applyDirichletValues(m_bc, 0, disc_vec);
  disc_vec->updateDependentDirichletValues();

  disc_vec->getDisc()->getMesh()->writeVtkFiles("dirichlet_test");

  

  for (int vol_block=0; vol_block < 2; ++vol_block)
  {
    auto& dofs = dof_numbering->getDofs(vol_block);
    auto& arr = disc_vec->getArray(vol_block);
    auto& coords = disc->getVolDisc(vol_block)->vol_group.coords;

    for (unsigned int i=0; i < dofs.shape()[0]; ++i)
      for (unsigned int j=0; j < dofs.shape()[1]; ++j)
      {
        auto& dofs = dof_numbering->getDofs(vol_block);

        //if (!dof_numbering->isDofActive(dofs[i][j]))
        //{
        //  std::cout << "non-active dof at " << coords[i][j][0] << ", " << coords[i][j][1] << ", " << coords[i][j][2] << std::endl;
        //  std::cout << "block " << vol_block << ", el " << i << ", node " << j << std::endl;
        //}
        Real val_ex = func(coords[i][j][0], coords[i][j][1], coords[i][j][2]);
        EXPECT_EQ(arr[i][j], val_ex);        
      }
  }
}

TEST_F(DiscVectorTester, AssignmentOperator)
{
  auto disc_vec2 = makeDiscVector(disc);
  disc_vec->setFunc(func);
  disc_vec2->set(0);

  disc_vec2->markArrayModified();

  *disc_vec2 = *disc_vec;

  EXPECT_TRUE(disc_vec2->isArrayCurrent());
  EXPECT_TRUE(disc_vec2->isVectorCurrent());

  auto& vec = disc_vec->getVector();
  auto& vec2 = disc_vec2->getVector();
  for (int i=0; i < disc_vec->getNumDofs(); ++i)
    EXPECT_EQ(vec[i], vec2[i]);

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto& arr = disc_vec->getArray(i);
    auto& arr2 = disc_vec2->getArray(i);
    for (int el=0; el < arr.shape()[0]; ++el)
      for (int j=0; j < arr.shape()[1]; ++j)
        EXPECT_EQ(arr[el][j], arr2[el][j]);
  }
}