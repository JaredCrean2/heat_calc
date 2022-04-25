#include "mesh/apfMDSField.h"

#ifdef MESH_USE_MDS_NUMBERING

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

}

#endif