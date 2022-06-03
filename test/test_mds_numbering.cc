#include "gtest/gtest.h"
#include "mesh/mesh_helper.h"
#include "test_helper.h"
#include <apfMDS.h>
#include <apfNumbering.h>
#include "mesh/apfMDSField.h"
#include "mesh/apfMDSFieldReduction.h"
#include "mesh/mesh_generator.h"
#include "mesh/reference_element_interface.h"
#include "mesh/reference_element_apf.h"

#ifdef MESH_USE_MDS_NUMBERING

namespace {
class MdsNumberingTester : public testing::Test
{
  protected:
    ~MdsNumberingTester()
    {
      if (numbering)
        apf::destroyNumbering(numbering);

      if (mesh)
        apf::destroyMesh(mesh);
    }

    void setup()
    {
      Mesh::MeshSpec spec;
      spec.nx = 5; spec.ny = 5; spec.nz = 5;
      setup(spec);
    }

    void setup(const Mesh::MeshSpec& mesh_spec)
    {
      if (numbering)
        apf::destroyNumbering(numbering);

      if (mesh)
        apf::destroyMesh(mesh);

      int comm_size;
      MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
      mesh = Mesh::make_parallel_mesh(mesh_spec, comm_size, &(Mesh::identity));

      apf::reorderMdsMesh(mesh);
  
      ref_el = reference_element::getLagrangeHexReferenceElement(3);
      fshape = apf::getHexFieldShape(ref_el);


      numbering = apf::createNumberingMDS(mesh, "foo", fshape.get(), ncomp);
    }

  apf::Mesh2* mesh = nullptr;
  reference_element::REPtr ref_el;
  std::shared_ptr<apf::FieldShapeRefEl> fshape;
  apf::ApfMDSNumbering* numbering = nullptr;
  int ncomp = 4;
};


class MdsFieldTester : public testing::Test
{
  protected:
    ~MdsFieldTester()
    {
      if (mesh)
        apf::destroyMesh(mesh);
    }

    void setup(bool add_ghost_layer=false)
    {
      Mesh::MeshSpec spec;
      spec.nx = 5; spec.ny = 5; spec.nz = 5;
      setup(spec, add_ghost_layer);
    }

    void setup(const Mesh::MeshSpec& mesh_spec, bool add_ghost_layer)
    {
      if (mesh)
        apf::destroyMesh(mesh);

      //auto generator = Mesh::make_mesh_generator(mesh_spec, &(Mesh::identity));
      //mesh = generator.generate();
      int comm_size;
      MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
      mesh = Mesh::make_parallel_mesh(mesh_spec, comm_size, &(Mesh::identity));
      if (add_ghost_layer)
        Mesh::createGhostLayer(mesh);

      apf::reorderMdsMesh(mesh);
  
      ref_el = reference_element::getLagrangeHexReferenceElement(2);
      field = std::make_shared<fast_field::ApfMDSField<double>>(mesh, "field_test", mesh->getShape(), ncomp);
    }

  apf::Mesh2* mesh = nullptr;
  reference_element::REPtr ref_el;
  std::shared_ptr<fast_field::ApfMDSField<double>> field = nullptr;
  int ncomp = 4;
};

}

TEST_F(MdsNumberingTester, IntializedUnnumbered)
{
  setup();
  for (int dim=0; dim < 3; ++dim)
  {
    int nnodes = fshape->countNodesOn(dim);
    if (nnodes == 0)
      continue;

    auto it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
      for (int j=0; j < nnodes; ++j)
        for (int comp=0; comp < ncomp; ++comp)
        {
          EXPECT_FALSE(apf::isNumbered(numbering, e, j, comp));
        }
    mesh->end(it);
  }
}


TEST_F(MdsNumberingTester, Number)
{
  setup();
  int idx = 0;
  for (int dim=0; dim < 3; ++dim)
  {
    int nnodes = fshape->countNodesOn(dim);
    if (nnodes == 0)
      continue;

    auto it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
      for (int j=0; j < nnodes; ++j)
        for (int comp=0; comp < ncomp; ++comp)
        {
          apf::number(numbering, e, j, comp, idx++);
        }
    mesh->end(it);
  }

  idx = 0;
  for (int dim=0; dim < 3; ++dim)
  {
    int nnodes = fshape->countNodesOn(dim);
    if (nnodes == 0)
      continue;

    auto it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
      for (int j=0; j < nnodes; ++j)
        for (int comp=0; comp < ncomp; ++comp)
        {
          EXPECT_EQ(apf::getNumber(numbering, e, j, comp), idx++);
        }
    mesh->end(it);
  }
}


TEST_F(MdsNumberingTester, Fixed)
{
  setup();
  for (int dim=0; dim < 3; ++dim)
  {
    int nnodes = fshape->countNodesOn(dim);
    if (nnodes == 0)
      continue;

    auto it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
      for (int j=0; j < nnodes; ++j)
        for (int comp=0; comp < ncomp; ++comp)
        {
          apf::fix(numbering, e, j, comp, true);
        }
    mesh->end(it);
  }

  for (int dim=0; dim < 3; ++dim)
  {
    int nnodes = fshape->countNodesOn(dim);
    if (nnodes == 0)
      continue;

    auto it = mesh->begin(dim);
    apf::MeshEntity* e;
    while ( (e = mesh->iterate(it)) )
      for (int j=0; j < nnodes; ++j)
        for (int comp=0; comp < ncomp; ++comp)
        {
          EXPECT_TRUE(apf::isFixed(numbering, e, j, comp));
          EXPECT_FALSE(apf::isNumbered(numbering, e, j, comp));
        }
    mesh->end(it);
  }
}

TEST_F(MdsNumberingTester, Accessors)
{  
  setup();
  EXPECT_EQ(apf::getShape(numbering), fshape.get());
  EXPECT_EQ(apf::getMesh(numbering), mesh);
  EXPECT_EQ(apf::countComponents(numbering), ncomp);
}

TEST_F(MdsNumberingTester, Synchronize)
{
  //setup(false);
  setup();
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  int valid_val   = std::numeric_limits<int>::max() - 2;
  int invalid_val = std::numeric_limits<int>::max() - 1;
  for (int dim=0; dim <= 3; ++dim)
  {
    if (!mesh->getShape()->hasNodesIn(dim))
      continue;

    it = mesh->begin(dim);
    while ( (e = mesh->iterate(it)) )
    {
      int etype = mesh->getType(e);
      for (int node=0; node < mesh->getShape()->countNodesOn(etype); ++node)
      {
        for (int component=0; component < ncomp; ++component)
          (*numbering)(e, node, component) = mesh->isOwned(e) ? valid_val : invalid_val;
      }
    }
    mesh->end(it);
  }

  apf::synchronize(numbering);

  for (int dim=0; dim <= 3; ++dim)
  {
    if (!mesh->getShape()->hasNodesIn(dim))
      continue;

    it = mesh->begin(dim);
    while ( (e = mesh->iterate(it)) )
    {
      int etype = mesh->getType(e);
      for (int node=0; node < mesh->getShape()->countNodesOn(etype); ++node)
      {
        for (int component=0; component < ncomp; ++component)
        {
          EXPECT_EQ((*numbering)(e, node, component), valid_val);
        }
      }
    }
    mesh->end(it);
  }
}


TEST_F(MdsFieldTester, Reduction)
{
  setup(true);
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  for (int dim=0; dim <= 3; ++dim)
  {
    if (!mesh->getShape()->hasNodesIn(dim))
      continue;

    it = mesh->begin(dim);
    while ( (e = mesh->iterate(it)) )
    {
      int etype = mesh->getType(e);
      for (int node=0; node < mesh->getShape()->countNodesOn(etype); ++node)
      {
        apf::Vector3 pt;
        mesh->getPoint(e, node, pt);
        double val = pt.x() + 2*pt.y() + 3*pt.z();
        for (int component=0; component < ncomp; ++component)
          (*field)(e, node, component) = val + component;
      }
    }
    mesh->end(it);
  }

  auto reducer = fast_field::make_field_reduction(*field, std::plus<double>());
  reducer.execute();

  for (int dim=0; dim <= 3; ++dim)
  {
    if (!mesh->getShape()->hasNodesIn(dim))
      continue;

    it = mesh->begin(dim);
    while ( (e = mesh->iterate(it)) )
    {
      int etype = mesh->getType(e);

      apf::Copies copies;
      if (mesh->isShared(e))
      {
        apf::Copies remotes;
        mesh->getRemotes(e, remotes);
        for (auto& p : remotes)
          copies.insert(p);
      }

      if (mesh->isGhost(e) || mesh->isGhosted(e))
      {
        apf::Copies ghosts;
        mesh->getGhosts(e, ghosts);
        for (auto& p : ghosts)
          copies.insert(p);
      }

      int count = copies.size() + 1;


      for (int node=0; node < mesh->getShape()->countNodesOn(etype); ++node)
      {
        apf::Vector3 pt;
        mesh->getPoint(e, node, pt);
        double val = pt.x() + 2*pt.y() + 3*pt.z();
        for (int component=0; component < ncomp; ++component)
        {
          EXPECT_EQ((*field)(e, node, component), count*(val + component));
        }
      }
    }
    mesh->end(it);
  }
}

#endif