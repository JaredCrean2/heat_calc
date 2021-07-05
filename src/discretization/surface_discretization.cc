#include "discretization/surface_discretization.h"

void computeNormals(const SurfaceDiscretization& disc, ArrayType<Real, 3> normals)
{
  using range = boost::multi_array_types::index_range;
  normals.resize(boost::extents[disc.getNumFaces()][disc.getNumSolPtsPerFace()][3]);

  for (int i=0; i < disc.getNumFaces(); ++i)
  {
    

  }



}
