#include <vector>
#include <cassert>
#include "ProjectDefs.h"

// class for evaluating Lagrange polynomials and their derivatives
class LagrangeBasis
{
  public: 
    using Index = std::vector<double>::size_type;

    LagrangeBasis(const std::vector<double>& pts);
    // evaluate the pth polynomial at point x
    double evalPoly(Index p, double x);

    double evalPolyDeriv(Index p, double x);

    Index getNumPoints() const { return m_pts.size();}

  private:
    std::vector<double> m_pts;
    std::vector<double> m_denoms;

};


// class that caches Lagrange polynomial values and derivatives for
// 3D tensor-product point sets
class LagrangeEvaluatorTP
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTP(const std::vector<double> pts_in,
                        const std::vector<double> pts_out) :
      m_vals(boost::extents[pts_out.size()][pts_in.size()]), 
      m_derivs(boost::extents[pts_out.size()][pts_in.size()])
    {
      LagrangeBasis basis(pts_in);
      for (Index i=0; i < pts_out.size(); ++i)
        for (Index j=0; j < pts_in.size(); ++j)
        {
          m_vals[i][j] = basis.evalPoly(j, pts_out[i]);
          m_derivs[i][j] = basis.evalPolyDeriv(j, pts_out[i]);
        }
    }

    // interpolates the values (defined at pts_int) to pts_out
    // It is recommended that T have the BOOST_RESTRICT qualifier for best performance
    template<typename Array>
    void interpolateVals(Array& vals_in, Array& vals_out)
    {
      assert(vals_in.num_dimensions()   == 3);
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
    void interpolateDerivs(Array3& vals_in, Array4& vals_out)
    {
      assert(vals_in.num_dimensions()   == 3);
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

    // TODO: check ordering of in vs out
    double getInterpolantValue(const Index i_in,  const Index j_in,  const Index k_in,
                               const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    double getInterpolantDx(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_derivs[i_out][i_in] * m_vals[j_out][j_in] * m_vals[k_out][k_in];
    }

    double getInterpolantDy(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_derivs[j_out][j_in] * m_vals[k_out][k_in];
    }

    double getInterpolantDz(const Index i_in,  const Index j_in,  const Index k_in,
                            const Index i_out, const Index j_out, const Index k_out)
    {
      return m_vals[i_out][i_in] * m_vals[j_out][j_in] * m_derivs[k_out][k_in];
    }

    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

  private:
    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    ArrayType<double, 2> m_vals;
    ArrayType<double, 2> m_derivs;
};


// non-tensor product output version
class LagrangeEvaluatorTPIn
{
  public:
    using Index = LagrangeBasis::Index;

    template <typename Array2D>
    LagrangeEvaluatorTPIn(const std::vector<Real>& pts_in,
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
    void interpolateVals(const Array3D& vals_in, Array1D& vals_out)
    {
      assert(vals_in.num_dimensions()   == 3);
      assert(vals_out.num_dimensions() == 1);
      for (int i=0; i < 3; ++i)
        assert(vals_in.shape()[i] == getNumPointsIn());

      assert(vals_out.shape()[0] == getNumPointsOut());

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
    void interpolateDerivs(const Array3D& vals_in, Array2D& vals_out)
    {
      assert(vals_in.num_dimensions()   == 3);
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



