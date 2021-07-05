#include "mesh/mesh.h"
#include <memory>

class SurfaceDiscretization
{
  public:
    SurfaceDiscretization(const FaceGroup& face_group, const std::vector<std::shared_ptr<VolumeGroup>> volume_groups) :
      face_group(face_group),
      volume_groups(volume_groups)
    {
      //TODO: set quadrature domain to be correct
      computeNormals(*this, normals);
    }

    int getNumFaces() const {return face_group.getNumFaces();}

    int getNumSolPtsPerFace() const { return face_group.getNumSolPtsPerFace()}

    ArrayType<Real, 3> normals;
    const FaceGroup& face_group;
    Quadrature quad;
    std::vector<std::shared_ptr<VolumeGroup>> volume_groups;
};


template <typename Array2D, typename ArrayNormalsXi, typename ArrayNormalsX>
void computeNormalNode(const ArrayIn& dxidx,
                       const ArrayNormalsXi& normals_xi,
                       ArrayNormalsX& normals_x)
{
  for (int i=0; i < 3; ++i)
    normals_x[i] = dxidx[0][i] * normals_xi[0] + 
                   dxidx[1][i] * normals_xi[1] +
                   dxidx[2][i] * normals_xi[2];
}
