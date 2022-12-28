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
    PostProcessorSurfaceIntegralAverage(const std::vector<SurfDiscPtr>& surfs, const std::string& name, Tfunc func) :
      m_surfs(surfs),
      m_name(name),
      m_func(func)
    {}

    int numValues() const override { return 1; }

    std::vector<std::string> getNames() const override { return {m_name}; }

    std::vector<double> getValues(DiscVectorPtr u, double t) override
    {
      if (!u->isArrayCurrent())
        u->syncArrayToVector();
       
      Real total_val = 0, total_area = 0;
      for (auto& surf : m_surfs)
      {

        ArrayType<Real, 1> vals(boost::extents[surf->getNumQuadPtsPerFace()]);

        for (int i=0; i < surf->getNumFaces(); ++i)
        {
          const Mesh::FaceSpec& spec = surf->face_group.faces[i];
          auto& u_arr = u->getArray(spec.vol_group);
          auto u_i = u_arr[boost::indices[spec.el_group][range()]];
          surf->interp_vsq_flat[spec.face].interpolateVals(u_i, vals);

          for (int j=0; j < surf->getNumQuadPtsPerFace(); ++j)
            vals[j] = m_func(vals[j]);

          total_val += integrateFaceScalar(surf, i, vals);

          for (int j=0; j < surf->getNumQuadPtsPerFace(); ++j)
            vals[j] = 1;

          total_area += integrateFaceScalar(surf, i, vals);
        }
      }

      return {total_val / total_area};
    }


  protected:
    std::vector<SurfDiscPtr> m_surfs;
    std::string m_name;
    Tfunc m_func;
};

template <typename T>
PostProcessorPtr makePostProcessorSurfaceIntegralAverage(SurfDiscPtr surf, const std::string& name, T func)
{
  return std::make_shared<PostProcessorSurfaceIntegralAverage<T>>(std::vector<SurfDiscPtr>{surf}, name, func);
}

template <typename T>
PostProcessorPtr makePostProcessorSurfaceIntegralAverage(const std::vector<SurfDiscPtr>& surfs, const std::string& name, T func)
{
  return std::make_shared<PostProcessorSurfaceIntegralAverage<T>>(surfs, name, func);
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