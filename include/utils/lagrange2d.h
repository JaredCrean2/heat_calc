#ifndef LAGRANGE2D_H
#define LAGRANGE2D_H

#include "utils/lagrange.h"

// class that caches Lagrange polynomial values and derivatives for
// 2D tensor-product point sets
class LagrangeEvaluatorTP2D
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTP2D(const std::vector<double> pts_in,
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
    void interpolateVals(Array& vals_in, Array& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 2);
      assert(vals_out.num_dimensions() == 2);
      for (int i=0; i < 2; ++i)
      {
        assert(vals_in.shape()[i] == getNumPointsIn());
        assert(vals_out.shape()[i] == getNumPointsOut());
      }


      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          {
            vals_out[i_out][j_out] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                  vals_out[i_out][j_out] += 
                    m_vals[i_out][i_in]*m_vals[j_out][j_in]*vals_in[i_in][j_in];
          }
    }

    // computes the derivatives of the values (defined at pts_in) at pts_out
    // vals_out is npts_out x npts_out x npts_out x 3
    template<typename Array2, typename Array3>
    void interpolateDerivs(Array2& vals_in, Array3& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 2);
      assert(vals_out.num_dimensions()  == 3);
      assert(vals_out.shape()[2]        == 2);
      for (int i=0; i < 2; ++i)
      {
        assert(vals_in.shape()[i] == getNumPointsIn());
        assert(vals_out.shape()[i] == getNumPointsOut());
      }


      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
        {
            vals_out[i_out][j_out][0] = 0;
            vals_out[i_out][j_out][1] = 0;
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
              {
                vals_out[i_out][j_out][0] += 
                  m_derivs[i_out][i_in]*m_vals[j_out][j_in]*vals_in[i_in][j_in];
                vals_out[i_out][j_out][1] += 
                  m_vals[i_out][i_in]*m_derivs[j_out][j_in]*vals_in[i_in][j_in];
                vals_out[i_out][j_out][2] += 
                  m_vals[i_out][i_in]*m_vals[j_out][j_in]*vals_in[i_in][j_in];

              }
        }
    }

    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

  private:
    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    ArrayType<double, 2> m_vals;
    ArrayType<double, 2> m_derivs;
};


// class that interpolates data that has tensor product structure, but 
// is given in a flat (1D) format
class LagrangeEvaluatorTPFlat2D
{
  public:
    using Index = LagrangeBasis::Index;

    LagrangeEvaluatorTPFlat2D(const std::vector<double>& pts_in,
                        const std::vector<double>& pts_out,
                        const ArrayType<LocalIndex, 2>& nodemap_in,
                        const ArrayType<LocalIndex, 2>& nodemap_out) :
      m_vals(boost::extents[pts_out.size()][pts_in.size()]), 
      m_derivs(boost::extents[pts_out.size()][pts_in.size()]),
      m_nodemap_in(nodemap_in),
      m_nodemap_out(nodemap_out)
    {
      assert(pts_in.size() == nodemap_in.shape()[0]);
      assert(pts_in.size() == nodemap_in.shape()[1]);
      assert(pts_out.size() == nodemap_out.shape()[0]);
      assert(pts_out.size() == nodemap_out.shape()[1]);

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
    void interpolateVals(Array& vals_in, Array& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 1);
      assert(vals_out.num_dimensions() == 1);
      assert(vals_in.shape()[0] == std::pow(getNumPointsIn(), 2));
      assert(vals_out.shape()[0] == std::pow(getNumPointsOut(), 2));


      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
          {
            vals_out[i_out][j_out] = 0;
            auto node_out = m_nodemap_out[i_out][j_out];
            for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
              for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
                  vals_out[node_out] += m_vals[i_out][i_in]*
                                 m_vals[j_out][j_in]*
                                 vals_in[m_nodemap_in[i_in][j_in]];
          }
    }

    // computes the derivatives of the values (defined at pts_in) at pts_out
    // vals_out is npts_out x npts_out x npts_out x 3
    template<typename Array1, typename Array2>
    void interpolateDerivs(Array1& vals_in, Array2& vals_out) const
    {
      assert(vals_in.num_dimensions()   == 1);
      assert(vals_in.shape()[0] == std::pow(getNumPointsIn(), 2));
      assert(vals_out.num_dimensions() == 2);
      assert(vals_out.shape()[0] == std::pow(getNumPointsOut(), 2));
      assert(vals_out.shape()[1]       == 3);

      for (Index i_out=0; i_out < getNumPointsOut(); ++i_out)
        for (Index j_out=0; j_out < getNumPointsOut(); ++j_out)
        {
          auto node_out = m_nodemap_out[i_out][j_out];
          vals_out[node_out][0] = 0;
          vals_out[node_out][1] = 0;
          for (Index i_in=0; i_in < getNumPointsIn(); ++i_in)
            for (Index j_in=0; j_in < getNumPointsIn(); ++j_in)
            {
              auto val = vals_in[m_nodemap_in[i_in][j_in]];
              vals_out[node_out][0] += 
                m_derivs[i_out][i_in]*m_vals[j_out][j_in]*val;
              vals_out[node_out][1] += 
                m_vals[i_out][i_in]*m_derivs[j_out][j_in]*val;
            }
        }
    }

    Index getNumPointsIn() const {return m_vals.shape()[1];}

    Index getNumPointsOut() const {return m_vals.shape()[0];}

  private:
    // TODO use BOOST_RESTRICT
    // basis function values and derivatives, npts_out x npts_in
    ArrayType<Real, 2> m_vals;
    ArrayType<Real, 2> m_derivs;
    const ArrayType<LocalIndex, 2>& m_nodemap_in;
    const ArrayType<LocalIndex, 2>& m_nodemap_out;
};



#endif
