#ifndef PHYSICS_POST_PROCESSOR
#define PHYSICS_POST_PROCESSOR

#include "post_processor_base.h"
#include "discretization/surface_discretization.h"
#include "discretization/NeumannBC.h"


namespace physics {

template <typename Tfunc>
class PostProcessorSurfaceIntegralAverage : public PostProcessorBase
{
  public: 
    PostProcessorSurfaceIntegralAverage(SurfDiscPtr surf, const std::string& name, Tfunc func) :
      m_surf(surf),
      m_name(name),
      m_func(func)
    {}

    int numValues() const override { return 1; }

    std::vector<std::string> getNames() const override { return {m_name}; }

    std::vector<double> getValues(DiscVectorPtr u, double t) override
    {
      if (!u->isArrayCurrent())
        u->syncArrayToVector();

      ArrayType<Real, 1> vals(boost::extents[m_surf->getNumQuadPtsPerFace()]);

      Real total_val = 0, total_area = 0;
      for (int i=0; i < m_surf->getNumFaces(); ++i)
      {
        const Mesh::FaceSpec& spec = m_surf->face_group.faces[i];
        auto& u_arr = u->getArray(spec.vol_group);
        auto u_i = u_arr[boost::indices[spec.el_group][range()]];
        m_surf->interp_vsq_flat[spec.face].interpolateVals(u_i, vals);

        for (int j=0; j < m_surf->getNumQuadPtsPerFace(); ++j)
          vals[j] = m_func(vals[j]);

        total_val += integrateFaceScalar(m_surf, i, vals);

        for (int j=0; j < m_surf->getNumQuadPtsPerFace(); ++j)
          vals[j] = 1;

        total_area += integrateFaceScalar(m_surf, i, vals);
      }

      return {total_val / total_area};
    }


  protected:
    SurfDiscPtr m_surf;
    std::string m_name;
    Tfunc m_func;
};

template <typename T>
PostProcessorPtr makePostProcessorSurfaceIntegralAverage(SurfDiscPtr surf, const std::string& name, T func)
{
  return std::make_shared<PostProcessorSurfaceIntegralAverage<T>>(surf, name, func);
}


class PostProcessorBCFlux : public PostProcessorBase
{
  public:

    PostProcessorBCFlux(const std::string& name, NeumannBCPtr bc) :
      m_name(name),
      m_bc(bc)
    {}

    int numValues() const override { return 1; }

    std::vector<std::string> getNames() const override { return {m_name}; }

    std::vector<double> getValues(DiscVectorPtr u, double t) override;

  private:
    std::string m_name;
    NeumannBCPtr m_bc;
};

}  // namespace

#endif