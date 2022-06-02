#include "mesh/apfMDSField.h"
#include "mesh/apfMDSFieldSynchronize.h"
#include "mesh/apfMDSFieldReduction.h"

#ifdef MESH_USE_MDS_NUMBERING

namespace fast_field {

// Setter must be callable as setter(MeshEntity*, int node, T* vals)
template <typename T, typename MeshWrapper, typename Setter>
void copyToApf(ApfMDSField<T, MeshWrapper>* field_in, Setter setter)
{
  apf::Mesh* m = field_in->getMesh();
  apf::FieldShape* fshape = field_in->getFieldShape();

  std::vector<T> vals(field_in->getNumComponents());
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  for (int dim=0; dim < 3; ++dim)
  {
    if (!fshape->hasNodesIn(dim))
      continue;

    it = m->begin(dim);
    while ( (e = m->iterate(it)) )
      for (int node=0; node < fshape->countNodesOn(m->getType(e)); ++node)
      {
        for (int c=0; c < field_in->getNumComponents(); ++c)
          vals[c] = (*field_in)(e, node, c);

        setter(e, node, vals.data());
      }
  }

}

} // namespace


namespace apf {

ApfMDSNumbering* createNumberingMDS(Mesh2* mesh, const char* name, FieldShape* shape, int component)
{
  return new ApfMDSNumbering(mesh, name, shape, component);
}

void destroyNumbering(ApfMDSNumbering* n)
{
  delete n;
}


FieldShape* getShape(ApfMDSNumbering* n)
{
  return n->getFieldShape();
}

const char* getName(ApfMDSNumbering* n)
{
  return n->getName().c_str();
}

Mesh* getMesh(ApfMDSNumbering* n)
{
  return n->getMesh();
}


void synchronize(ApfMDSNumbering* n)
{
  auto syncer = fast_field::make_field_synchronizer(*n);
  //auto syncer = fast_field::make_field_reduction(*n, fast_field::assign<int>());
  syncer.synchronize();
}

void copyToApf(ApfMDSNumbering* n_in, apf::Numbering* n_out)
{
  assertAlways(n_in->getFieldShape() == apf::getShape(n_out), "FieldShapes must be the same");

  auto setter = [&](apf::MeshEntity* e, int node, int* vals)
  {
    for (int c=0; c < n_in->getNumComponents(); ++c)
      apf::number(n_out, e, node, c, vals[c]);
  };

  copyToApf(n_in, setter);
}

}

#endif