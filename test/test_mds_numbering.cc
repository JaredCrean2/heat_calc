#include "gtest/gtest.h"
#include "test_helper.h"
#include <apfMDS.h>
#include <apfNumbering.h>
#include "mesh/apfMDSField.h"
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

      auto generator = Mesh::make_mesh_generator(mesh_spec, &(Mesh::identity));
      mesh = generator.generate();
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

}

TEST_F(MdsNumberingTester, IntializedUnnumbered)
{
  SERIAL_ONLY();

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
  SERIAL_ONLY();

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
  SERIAL_ONLY();

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
  SERIAL_ONLY();
  
  setup();
  EXPECT_EQ(apf::getShape(numbering), fshape.get());
  EXPECT_EQ(apf::getMesh(numbering), mesh);
  EXPECT_EQ(apf::countComponents(numbering), ncomp);
}

#endif