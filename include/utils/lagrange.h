#ifndef LAGRANGE_H
#define LAGRANGE_H

#include <vector>
#include <cassert>
#include <iostream>

#include "ProjectDefs.h"
#include "lagrange_basis.h"


// class that caches Lagrange polynomial values and derivatives for
// 3D tensor-product point sets
class LagrangeEvaluatorTPToTP
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTPToTP(const std::vector<Real> pts_in,
                            const std::vector<Real> pts_out) :
      m_vals(lagrange_memoizer.getValues(pts_in, pts_out)),
      m_derivs(lagrange_memoizer.getDerivs(pts_in, pts_out))
    {}

    // interpolates the values (defined at pts_int) to pts_out
    // It is recommended that T have the BOOST_RESTRICT qualifier for best performance
    template<typename ArrayIn, typename ArrayOut>
    void interpolateVals(ArrayIn& vals_in, ArrayOut& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 3);
      for (int i=0; i < 3; ++i)
      {
        assert(vals_in.shape()[i] == getNumPointsIn());
        assert(vals_out.shape()[i] == getNumPointsOut());
      }

      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumPointsOut(); ++k_out)
          {
            vals_out[i_out][j_out][k_out] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                //TODO: precompute first two lagrange polynomials here
                for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
                  vals_out[i_out][j_out][k_out] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*m_vals[k_out][k_in]*vals_in[i_in][j_in][k_in];
          }
    }

    // computes the derivatives of the values (defined at pts_in) at pts_out
    // vals_out is npts_out x npts_out x npts_out x 3
    template<typename Array3, typename Array4>
    void interpolateDerivs(Array3& vals_in, Array4& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 4);
      assert(vals_out.shape()[3]       == 3);
      for (int i=0; i < 3; ++i)
      {
        assert(vals_in.shape()[i] == getNumPointsIn());
        assert(vals_out.shape()[i] == getNumPointsOut());
      }

      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumPointsOut(); ++k_out)
          {
            vals_out[i_out][j_out][k_out][0] = 0;
            vals_out[i_out][j_out][k_out][1] = 0;
            vals_out[i_out][j_out][k_out][2] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                // TODO: precompute first two multiplications
                for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
                {
                  vals_out[i_out][j_out][k_out][0] += 
                    m_derivs[i_out][i_in]*m_vals[j_out][j_in]*m_vals[k_out][k_in]*vals_in[i_in][j_in][k_in];
                  vals_out[i_out][j_out][k_out][1] += 
                    m_vals[i_out][i_in]*m_derivs[j_out][j_in]*m_vals[k_out][k_in]*vals_in[i_in][j_in][k_in];
                  vals_out[i_out][j_out][k_out][2] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*m_derivs[k_out][k_in]*vals_in[i_in][j_in][k_in];

                }
          }
    }
/*
    // TODO: check ordering of in vs out
    Real getInterpolantValue(const Index i_in,  const Index j_in,  const Index k_in,
                               const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    Real getInterpolantDx(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_derivs[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    Real getInterpolantDy(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_derivs[j_out][j_in] * m_vals[k_out][k_in];
    }

    Real getInterpolantDz(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_derivs[k_out][k_in];
    }
*/
    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

  private:
    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    LagrangeMemoizer::RetType m_vals;
    LagrangeMemoizer::RetType m_derivs;
};



class LagrangeEvaluatorTPFlatToTP
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTPFlatToTP(const std::vector<Real> pts_in,
                                const std::vector<Real> pts_out,
                                const ArrayType<LocalIndex, 3>& nodemap_in) :
      m_vals(lagrange_memoizer.getValues(pts_in, pts_out)), 
      m_derivs(lagrange_memoizer.getDerivs(pts_in, pts_out)),
      m_nodemap_in(nodemap_in)
    {
      for (int i=0; i < 3; ++i)
        assert(pts_in.size() == nodemap_in.shape()[i]);
    }

    // interpolates the values (defined at pts_int) to pts_out
    // It is recommended that T have the BOOST_RESTRICT qualifier for best performance
    template<typename ArrayIn, typename ArrayOut>
    void interpolateVals(ArrayIn& vals_in, ArrayOut& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 1);
      assert(vals_out.num_dimensions() == 3);
      assert(vals_in.shape()[0] == getNumPointsIn());
      for (int i=0; i < 3; ++i)
        assert(vals_out.shape()[i] == getNumPointsOut());


      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumPointsOut(); ++k_out)
          {
            vals_out[i_out][j_out][k_out] = 0;
            for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
                //TODO: precompute first two lagrange polynomials here
                for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
                  vals_out[i_out][j_out][k_out] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*
                    m_vals[k_out][k_in]*
                    vals_in[m_nodemap_in[i_in][j_in][k_in]];
          }
    }

    // computes the derivatives of the values (defined at pts_in) at pts_out
    // vals_out is npts_out x npts_out x npts_out x 3
    template<typename Array1, typename Array4>
    void interpolateDerivs(Array1& vals_in, Array4& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 1);
      assert(vals_out.num_dimensions() == 4);
      assert(vals_in.shape()[0]        == getNumPointsIn());
      assert(vals_out.shape()[3]       == 3);
      for (int i=0; i < 3; ++i)
        assert(vals_out.shape()[i] == getNumPointsOut());


      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumPointsOut(); ++k_out)
          {
            vals_out[i_out][j_out][k_out][0] = 0;
            vals_out[i_out][j_out][k_out][1] = 0;
            vals_out[i_out][j_out][k_out][2] = 0;
            for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
                // TODO: precompute first two multiplications
                for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
                {
                  auto node = m_nodemap_in[i_in][j_in][k_in];
                  vals_out[i_out][j_out][k_out][0] += 
                    m_derivs[i_out][i_in]*m_vals[j_out][j_in]*m_vals[k_out][k_in]*vals_in[node];
                  vals_out[i_out][j_out][k_out][1] += 
                    m_vals[i_out][i_in]*m_derivs[j_out][j_in]*m_vals[k_out][k_in]*vals_in[node];
                  vals_out[i_out][j_out][k_out][2] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*m_derivs[k_out][k_in]*vals_in[node];

                }
          }
    }

    Index getNumPointsIn() const {return m_nodemap_in.num_elements();}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

  private:

    Index getNumTPPointsIn() const {return m_vals.shape()[1];}

    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    LagrangeMemoizer::RetType m_vals;
    LagrangeMemoizer::RetType m_derivs;
    const ArrayType<LocalIndex, 3>& m_nodemap_in;
};



class LagrangeEvaluatorTPToTPFlat
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTPToTPFlat(const std::vector<Real> pts_in,
                                const std::vector<Real> pts_out,
                                const ArrayType<LocalIndex, 3>& nodemap_out
                                ) :
      m_vals(lagrange_memoizer.getValues(pts_in, pts_out)), 
      m_derivs(lagrange_memoizer.getDerivs(pts_in, pts_out)),
      m_nodemap_out(nodemap_out)
    {
      for (int i=0; i < 3; ++i)
        assert(pts_out.size() == nodemap_out.shape()[i]);
    }

    // interpolates the values (defined at pts_in) to pts_out
    // It is recommended that T have the BOOST_RESTRICT qualifier for best performance
    template<typename ArrayIn, typename ArrayOut>
    void interpolateVals(ArrayIn& vals_in, ArrayOut& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 1);
      assert(vals_out.shape()[0]       == getNumPointsOut());
      for (int i=0; i < 3; ++i)
        assert(vals_in.shape()[i] == getNumPointsIn());


      for (Index i_out=0; i_out < getNumTPPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumTPPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumTPPointsOut(); ++k_out)
          {
            auto node = m_nodemap_out[i_out][j_out][k_out];
            vals_out[node] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                //TODO: precompute first two lagrange polynomials here
                for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
                  vals_out[node] += m_vals[i_out][i_in]*m_vals[j_out][j_in]*
                    m_vals[k_out][k_in]*vals_in[i_in][j_in][k_in];
          }
    }

    // computes the derivatives of the values
    template<typename Array3, typename Array2>
    void interpolateDerivs(Array3& vals_in, Array2& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 2);
      assert(vals_out.shape()[0]       == getNumPointsOut());
      assert(vals_out.shape()[1]       == 3);
      for (int i=0; i < 3; ++i)
        assert(vals_in.shape()[i] == getNumPointsIn());

      for (Index i_out=0; i_out < getNumTPPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumTPPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumTPPointsOut(); ++k_out)
          {
            auto node = m_nodemap_out[i_out][j_out][k_out];
            vals_out[node][0] = 0;
            vals_out[node][1] = 0;
            vals_out[node][2] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                // TODO: precompute first two multiplications
                for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
                {
                  vals_out[node][0] += m_derivs[i_out][i_in]*
                    m_vals[j_out][j_in]*m_vals[k_out][k_in]*
                    vals_in[i_in][j_in][k_in];
                  vals_out[node][1] += m_vals[i_out][i_in]*
                    m_derivs[j_out][j_in]*m_vals[k_out][k_in]*
                    vals_in[i_in][j_in][k_in];
                  vals_out[node][2] += m_vals[i_out][i_in]*
                    m_vals[j_out][j_in]*m_derivs[k_out][k_in]*
                    vals_in[i_in][j_in][k_in];
                }
          }
    }

    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_nodemap_out.num_elements();}

  private:

    Index getNumTPPointsOut() const {return m_vals.shape()[0];}

    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    LagrangeMemoizer::RetType m_vals;
    LagrangeMemoizer::RetType m_derivs;
    const ArrayType<LocalIndex, 3>& m_nodemap_out;
};


class LagrangeEvaluatorTPFlatToTPFlat
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTPFlatToTPFlat(const std::vector<Real>& pts_in,
                        const std::vector<Real>& pts_out,
                        const ArrayType<LocalIndex, 3>& nodemap_in,
                        const ArrayType<LocalIndex, 3>& nodemap_out) :
      m_vals(lagrange_memoizer.getValues(pts_in, pts_out)), 
      m_derivs(lagrange_memoizer.getDerivs(pts_in, pts_out)),
      m_nodemap_in(nodemap_in),
      m_nodemap_out(nodemap_out)
    {
      for (int i=0; i < 3; ++i)
      {
        assert(pts_in.size() == nodemap_in.shape()[i]);
        assert(pts_out.size() == nodemap_out.shape()[i]);
      }
    }

    // interpolates the values (defined at pts_int) to pts_out
    // It is recommended that T have the BOOST_RESTRICT qualifier for best performance
    template<typename ArrayIn, typename ArrayOut>
    void interpolateVals(ArrayIn& vals_in, ArrayOut& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 1);
      assert(vals_out.num_dimensions() == 1);
      assert(vals_in.shape()[0]        == getNumPointsIn());
      assert(vals_out.shape()[0]       == getNumPointsOut());


      for (Index i_out=0; i_out < getNumTPPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumTPPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumTPPointsOut(); ++k_out)
          {
            auto node_out = m_nodemap_out[i_out][j_out][k_out];
            vals_out[node_out] = 0;
            for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
                //TODO: precompute first two lagrange polynomials here
                for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
                  vals_out[node_out] += m_vals[i_out][i_in]*
                                 m_vals[j_out][j_in]*m_vals[k_out][k_in]*
                                 vals_in[m_nodemap_in[i_in][j_in][k_in]];
          }
    }

    // computes the derivatives of the values
    template<typename Array1, typename Array2>
    void interpolateDerivs(Array1& vals_in, Array2& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 1);
      assert(vals_out.num_dimensions() == 2);
      assert(vals_in.shape()[0]        == getNumPointsIn());
      assert(vals_out.shape()[0]       == getNumPointsOut());
      assert(vals_out.shape()[1]       == 3);


      for (Index i_out=0; i_out < getNumTPPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumTPPointsOut(); ++j_out)
          for (Index k_out=0; k_out < getNumTPPointsOut(); ++k_out)
          {
            auto node_out = m_nodemap_out[i_out][j_out][k_out];
            vals_out[node_out][0] = 0;
            vals_out[node_out][1] = 0;
            vals_out[node_out][2] = 0;
            for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
                // TODO: precompute first two multiplications
                for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
                {
                  auto val = vals_in[m_nodemap_in[i_in][j_in][k_in]];
                  vals_out[node_out][0] += 
                    m_derivs[i_out][i_in]*m_vals[j_out][j_in]*m_vals[k_out][k_in]*val;
                  vals_out[node_out][1] += 
                    m_vals[i_out][i_in]*m_derivs[j_out][j_in]*m_vals[k_out][k_in]*val;
                  vals_out[node_out][2] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*m_derivs[k_out][k_in]*val;

                }
          }
    }

    Index getNumPointsIn() const {return m_nodemap_in.num_elements();}

    Index getNumPointsOut() const {return m_nodemap_out.num_elements();}

  private:
    Index getNumTPPointsIn() const {return m_vals.shape()[1];}

    Index getNumTPPointsOut() const {return m_vals.shape()[0];}

    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    LagrangeMemoizer::RetType m_vals;
    LagrangeMemoizer::RetType m_derivs;
    const ArrayType<LocalIndex, 3>& m_nodemap_in;
    const ArrayType<LocalIndex, 3>& m_nodemap_out;
};





// non-tensor product output version
class LagrangeEvaluatorTPToNonTP
{
  public:
    using Index = LagrangeBasis::Index;

    template <typename Array2D>
    LagrangeEvaluatorTPToNonTP(const std::vector<Real>& pts_in,
                               const Array2D& pts_out) :
      m_vals(boost::extents[pts_out.size()][pts_in.size()][3]), 
      m_derivs(boost::extents[pts_out.shape()[0]][pts_in.size()][3])
    {
      assert(pts_out.num_dimensions() == 2);
      LagrangeBasis basis(pts_in);
      for (Index i=0; i < pts_out.size(); ++i)
        for (Index j=0; j < pts_in.size(); ++j)
          for (Index d=0; d < 3; ++d)
          {
            m_vals[i][j][d]   = basis.evalPoly(j, pts_out[i][d]);
            m_derivs[i][j][d] = basis.evalPolyDeriv(j, pts_out[i][d]);
          }
    }

    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

    template <typename Array3D, typename Array1D>
    void interpolateVals(const Array3D& vals_in, Array1D& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 1);
      assert(vals_out.shape()[0] == getNumPointsOut());
      for (int i=0; i < 3; ++i)
        assert(vals_in.shape()[i] == getNumPointsIn());


      for (Index i_out =0; i_out < getNumPointsOut(); ++i_out)
      {
        vals_out[i_out] = 0;
        for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
          for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
            for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
              // TODO: precompute first two multiplications
              vals_out[i_out] += m_vals[i_out][i_in][0] *
                                 m_vals[i_out][j_in][1] *
                                 m_vals[i_out][k_in][2] *
                                 vals_in[i_in][j_in][k_in];
      }
    }

    template <typename Array3D, typename Array2D>
    void interpolateDerivs(const Array3D& vals_in, Array2D& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 3);
      assert(vals_out.num_dimensions() == 2);
      for (int i=0; i < 3; ++i)
        assert(vals_in.shape()[i] == getNumPointsIn());

      assert(vals_out.shape()[0] == getNumPointsOut());
      assert(vals_out.shape()[1] == 3);

      for (Index i_out =0; i_out < getNumPointsOut(); ++i_out)
      {
        vals_out[i_out][0] = 0; vals_out[i_out][1] = 0;
        vals_out[i_out][2] = 0;
        for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
          for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
            for (Index k_in=0; k_in < getNumPointsIn(); ++k_in)
            {
              // TODO: precompute first two multiplications
              vals_out[i_out][0] += m_derivs[i_out][i_in][0] *
                                      m_vals[i_out][j_in][1] *
                                      m_vals[i_out][k_in][2] *
                                      vals_in[i_in][j_in][k_in];

              vals_out[i_out][1] +=   m_vals[i_out][i_in][0] *
                                    m_derivs[i_out][j_in][1] *
                                      m_vals[i_out][k_in][2] *
                                      vals_in[i_in][j_in][k_in];

              vals_out[i_out][2] +=   m_vals[i_out][i_in][0] *
                                      m_vals[i_out][j_in][1] *
                                    m_derivs[i_out][k_in][2] *
                                      vals_in[i_in][j_in][k_in];
            }
      }
    }


  private:
    ArrayType<Real, 3> m_vals;
    ArrayType<Real, 3> m_derivs;
};

// non-tensor product output version
class LagrangeEvaluatorTPFlatToNonTP
{
  public:
    using Index = LagrangeBasis::Index;

    template <typename Array2D>
    LagrangeEvaluatorTPFlatToNonTP(const std::vector<Real>& pts_in,
                               const Array2D& pts_out,
                               const ArrayType<LocalIndex, 3> nodemap_in) :
      m_vals(boost::extents[pts_out.size()][pts_in.size()][3]), 
      m_derivs(boost::extents[pts_out.shape()[0]][pts_in.size()][3]),
      m_nodemap_in(nodemap_in)
    {
      for (int i=0; i < 3; ++i)
        assert(pts_in.size() == nodemap_in.shape()[i]);

      assert(pts_out.num_dimensions() == 2);
      LagrangeBasis basis(pts_in);
      for (Index i=0; i < pts_out.size(); ++i)
        for (Index j=0; j < pts_in.size(); ++j)
          for (Index d=0; d < 3; ++d)
          {
            m_vals[i][j][d]   = basis.evalPoly(j, pts_out[i][d]);
            m_derivs[i][j][d] = basis.evalPolyDeriv(j, pts_out[i][d]);
          }
    }

    Index getNumPointsIn() const {return m_nodemap_in.num_elements();}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

    template <typename ArrayIn, typename ArrayOut>
    void interpolateVals(const ArrayIn& vals_in, ArrayOut& vals_out) const
    {
      assert(vals_in.num_dimensions()  == 1);
      assert(vals_out.num_dimensions() == 1);
      assert(vals_in.shape()[0] == getNumPointsIn());
      assert(vals_out.shape()[0] == getNumPointsOut());

      for (Index i_out =0; i_out < getNumPointsOut(); ++i_out)
      {
        vals_out[i_out] = 0;
        for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
          for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
            for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
              // TODO: precompute first two multiplications
              vals_out[i_out] += m_vals[i_out][i_in][0] *
                                 m_vals[i_out][j_in][1] *
                                 m_vals[i_out][k_in][2] *
                                 vals_in[m_nodemap_in[i_in][j_in][k_in]];
      }
    }

    template <typename Array1D, typename Array2D>
    void interpolateDerivs(const Array1D& vals_in, Array2D& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 1);
      assert(vals_out.num_dimensions() == 2);
      assert(vals_in.shape()[0] == getNumPointsIn());
      assert(vals_out.shape()[0] == getNumPointsOut());
      assert(vals_out.shape()[1] == 3);

      for (Index i_out =0; i_out < getNumPointsOut(); ++i_out)
      {
        vals_out[i_out][0] = 0; 
        vals_out[i_out][1] = 0;
        vals_out[i_out][2] = 0;
        for (Index i_in=0; i_in < getNumTPPointsIn(); ++i_in)
          for (Index j_in=0; j_in < getNumTPPointsIn(); ++j_in)
            for (Index k_in=0; k_in < getNumTPPointsIn(); ++k_in)
            {
              auto val = vals_in[m_nodemap_in[i_in][j_in][k_in]];
              // TODO: precompute first two multiplications
              vals_out[i_out][0] += m_derivs[i_out][i_in][0] *
                                      m_vals[i_out][j_in][1] *
                                      m_vals[i_out][k_in][2] *
                                      val;

              vals_out[i_out][1] +=   m_vals[i_out][i_in][0] *
                                    m_derivs[i_out][j_in][1] *
                                      m_vals[i_out][k_in][2] *
                                      val;

              vals_out[i_out][2] +=   m_vals[i_out][i_in][0] *
                                      m_vals[i_out][j_in][1] *
                                    m_derivs[i_out][k_in][2] *
                                      val;
            }
      }
    }


  private:
    Index getNumTPPointsIn() const {return m_vals.shape()[1];}

    ArrayType<Real, 3> m_vals;
    ArrayType<Real, 3> m_derivs;
    const ArrayType<LocalIndex, 3> m_nodemap_in;
};


#endif
