#include "gtest/gtest.h"
#include "mesh_helper.h"
#include <set>

namespace {
  class DofNumberingTester : public ::testing::Test,
                             public StandardDiscSetup
{
  protected:
    DofNumberingTester()
    {
      setup();
    }
};

}  // namespace


TEST_F(DofNumberingTester, Total)
{
  // count unique dofs
  auto dof_numbering = disc->getDofNumbering();

  int num_dofs_ex = (spec.nx + 1) * (spec.ny + 1) * (spec.nz + 1);
  int num_active_dofs_ex = (spec.nx - 1) * (spec.ny - 1) * (spec.nz - 1);

  std::set<Index> dof_nums_set;
  std::set<Index> active_dof_nums_set;
  auto vol_disc = disc->getVolDisc(0);
  auto& dofs = dof_numbering->getDofs(vol_disc);

  for (int i=0; i < vol_disc->getNumElems(); ++i)
    for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
    {
      Index dof = dofs[i][j];
      EXPECT_GE(dof, 0);
      EXPECT_LT(dof, num_dofs_ex);
      dof_nums_set.insert(dof);
      if (dof_numbering->isDofActive(dof))
        active_dof_nums_set.insert(dof);
    }

  EXPECT_EQ(dof_nums_set.size(), static_cast<size_t>(num_dofs_ex));
  EXPECT_EQ(active_dof_nums_set.size(), static_cast<size_t>(num_active_dofs_ex));
}


TEST_F(DofNumberingTester, Uniqueness)
{
  // test that nodes that have the same dof number also have the same coordinates
  auto dof_numbering = disc->getDofNumbering();
  auto vol_disc = disc->getVolDisc(0);
  auto& dofs = dof_numbering->getDofs(vol_disc);

  std::map<Index, std::vector<ElementNode>> dof_copies;
  for (int i=0; i < vol_disc->getNumElems(); ++i)
    for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      dof_copies[dofs[i][j]].emplace_back(i, j);

  for (auto& dofs_p : dof_copies)
    if (dofs_p.second.size() > 1)
    {
      auto& copies_p = dofs_p.second;
      auto& node_first = copies_p.front();
      auto& coords = vol_disc->vol_group.coords;
      auto coords_first = coords[boost::indices[node_first.element][node_first.node][range()]];

      for (unsigned int j=1; j < copies_p.size(); ++j)
      {
        auto node_j = copies_p[j];
        for (int d=0; d < 3; ++d)
          EXPECT_NEAR(coords_first[d], coords[node_j.element][node_j.node][d], 1e-13);
      }
    }
}


TEST_F(DofNumberingTester, DirichletDofsCount)
{

  auto dof_numbering = disc->getDofNumbering();
  auto vol_disc = disc->getVolDisc(0);

  int num_dirichlet_ex = 2*(spec.nx + 1)*(spec.ny + 1) + 2*(spec.ny + 1)*(spec.nz+1) +
                         2*(spec.nx + 1)*(spec.nz + 1) - 
                         4*(spec.nx - 1) - 4*(spec.ny - 1) - 4*(spec.nz - 1) -
                         2*8;

  const auto& dirichlet_dofs = dof_numbering->getDirichletNodes(vol_disc);

  EXPECT_EQ(dirichlet_dofs.size(), static_cast<size_t>(num_dirichlet_ex));
}


TEST_F(DofNumberingTester, DirichletDofsUniqueness)
{

  auto dof_numbering = disc->getDofNumbering();
  auto vol_disc = disc->getVolDisc(0);

  const auto& dofs = dof_numbering->getDofs(vol_disc);
  const auto& dirichlet_dofs = dof_numbering->getDirichletNodes(vol_disc);

  std::set<Index> dofs_seen;
  for (auto& node : dirichlet_dofs)
  {
    Index dof = dofs[node.element][node.node];
    EXPECT_TRUE(dofs_seen.count(dof) == 0);
    dofs_seen.insert(dof);
  }
}

TEST_F(DofNumberingTester, DirichletDofsCoords)
{
  auto dof_numbering = disc->getDofNumbering();
  auto vol_disc = disc->getVolDisc(0);

  const auto& dirichlet_dofs = dof_numbering->getDirichletNodes(vol_disc);

  for (auto& node : dirichlet_dofs)
  {
    auto coords = vol_disc->vol_group.coords[boost::indices[node.element][node.node][range()]];
    bool on_surface = false;
    for (int d=0; d < 3; ++d)
      on_surface |= std::abs(coords[d] - mesh_dim_mins[d]) < 1e-13 ||
                    std::abs(coords[d] - mesh_dim_maxs[d]) < 1e-13;
    EXPECT_TRUE(on_surface);
  }
}


TEST_F(DofNumberingTester, DirichletDofVolumeCoords)
{
  setup(5, 3);
  auto dof_numbering = disc->getDofNumbering();
  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    auto vol_disc = disc->getVolDisc(i);
    auto dofs_i = dof_numbering->getDofs(i);
    ArrayType<Real, 2> sol_coords(boost::extents[vol_disc->getNumSolPtsPerElement()][3]);
    for (int el=0; el < vol_disc->getNumElems(); ++el)
    {
      vol_disc->getVolumeSolCoords(el, sol_coords);
      for (int j=0; j < vol_disc->getNumSolPtsPerElement(); ++j)
      {
        auto dof = dofs_i[el][j];
        auto coords_j = sol_coords[boost::indices[j][range()]];

        bool on_surface = false;
        for (int d=0; d < 3; ++d)
          on_surface |= std::abs(coords_j[d] - mesh_dim_mins[d]) < 1e-13 ||
                        std::abs(coords_j[d] - mesh_dim_maxs[d]) < 1e-13;   

        //std::cout << "el, j = " << el << ", " << j << std::endl;
        //std::cout << std::boolalpha << "coords_j = " << coords_j[0] << ", " << coords_j[1] << ", " << coords_j[2] << ", is on surface = " << on_surface << std::endl;
        if (dof_numbering->isDofActive(dof))
          EXPECT_FALSE(on_surface);
        else
          EXPECT_TRUE(on_surface);     

      }
    }
  }

}