#ifndef MESH_REFERENCE_ELEMENT_H
#define MESH_REFERENCE_ELEMENT_H

#include "ProjectDefs.h"
#include "mesh/mesh_input.h"  // for LocalIndex etc.
#include "apfMesh.h"

namespace Mesh {

class ReferenceElement
{
  public:
    ReferenceElement(const int degree): degree(degree) {}
    virtual ~ReferenceElement() {}

    virtual const ArrayType<LocalIndex, 3>& getTensorProductMap() const = 0;

    virtual const std::vector<Real>& getTensorProductXi() const = 0;

    virtual const ArrayType<Real, 2>& getNodeXi() const = 0;

    virtual const ArrayType<Real, 2>& getNormals() const = 0;

    // given the xi coordinates on a face, computes the xi coordinates of
    // the same point in the elements xi coordinates
    virtual void computeElementXi(const int face, const Real* xi_face, Real* xi_element) const = 0;

    // get total number of nodes
    virtual int getNumNodes() const = 0;

    // get number of nodes in each direction
    virtual int getNumNodesTensorProduct() const = 0;

    int getDegree() const { return degree;}

    virtual int getNumFaces() const = 0;

    virtual std::pair<Real, Real> getXiRange() const { return std::pair<Real, Real>(0, 1);}

  protected:
    int degree;    
};

class ReferenceElementTensorProduct : public ReferenceElement
{
  public:
    ReferenceElementTensorProduct(const int degree) :
      ReferenceElement(degree)
    {}

    virtual ~ReferenceElementTensorProduct() {}

    const ArrayType<Real, 2>& getNodeXi() const override;

    int getNumNodes() const override;

    int getNumNodesTensorProduct() const override { return getTensorProductXi().size();}
};


class HexReferenceElement : public ReferenceElementTensorProduct
{
  public:
    HexReferenceElement(const int degree) :
      ReferenceElementTensorProduct(degree)
    {}

    const ArrayType<LocalIndex, 3>& getTensorProductMap() const override;

    const std::vector<Real>& getTensorProductXi() const override;

    const ArrayType<Real, 2>& getNormals() const override;

    // given the xi coordinates on a face, computes the xi coordinates of
    // the same point in the elements xi coordinates
    void computeElementXi(const int face, const Real* xi_face, Real* xi_element) const override;

    int getNumFaces() const override { return 6;}
};


ReferenceElement* getReferenceElement(apf::Mesh::Type type, const int degree);

}

#endif
