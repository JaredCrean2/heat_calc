#ifndef PHYSICS_POST_PROCESSOR
#define PHYSICS_POST_PROCESSOR

#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/bc_defs.h"
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

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override
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

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    std::string m_name;
    NeumannBCPtr m_bc;
};


class PostProcessorAirWindSkyBCFlux : public PostProcessorBCFlux
{
  public:

    PostProcessorAirWindSkyBCFlux(const std::string& name, std::shared_ptr<Heat::AirWindSkyNeumannBC> bc,
                                  Heat::HeatEquationSolar* heat_eqn_solar) :
      PostProcessorBCFlux(name, bc),
      m_heat_eqn_solar(heat_eqn_solar)
    {}

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    Heat::HeatEquationSolar* m_heat_eqn_solar;
};


class PostProcessorCombinedAirWindSkyBCFlux : public PostProcessorBase
{
  public:

    PostProcessorCombinedAirWindSkyBCFlux(const std::string& name_prefix, std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> bc,
                                  Heat::HeatEquationSolar* heat_eqn_solar) :
      m_name_prefix(name_prefix),
      m_bc(bc),
      m_heat_eqn_solar(heat_eqn_solar)
    {}

    // returns number of values this postprocessor returns
    int numValues() const override;

    std::vector<std::string> getNames() const override;

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    std::string m_name_prefix;
    std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> m_bc;
    Heat::HeatEquationSolar* m_heat_eqn_solar;    
};

}  // namespace

#endif