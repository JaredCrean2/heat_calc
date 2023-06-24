#ifndef PHYSICS_POST_PROCESSOR
#define PHYSICS_POST_PROCESSOR

#include "ProjectDefs.h"
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
    PostProcessorSurfaceIntegralAverage(const std::vector<SurfDiscPtr>& surfs, const std::string& name, Tfunc func, MPI_Comm comm) :
      m_surfs(surfs),
      m_name(name),
      m_func(func),
      m_comm(comm)
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

          total_val += surf->getFaceWeight(i) * integrateFaceScalar(surf, i, vals);

          for (int j=0; j < surf->getNumQuadPtsPerFace(); ++j)
            vals[j] = 1;

          total_area += surf->getFaceWeight(i) * integrateFaceScalar(surf, i, vals);
        }
      }

      std::array<Real, 2> vals_local = {total_val, total_area}, vals_global;
      MPI_Allreduce(vals_local.data(), vals_global.data(), 2, REAL_MPI_DATATYPE, MPI_SUM, m_comm);

      return {vals_global[0] / vals_global[1]};
    }


  protected:
    std::vector<SurfDiscPtr> m_surfs;
    std::string m_name;
    Tfunc m_func;
    MPI_Comm m_comm;
};

template <typename T>
PostProcessorPtr makePostProcessorSurfaceIntegralAverage(SurfDiscPtr surf, const std::string& name, T func, MPI_Comm comm)
{
  return std::make_shared<PostProcessorSurfaceIntegralAverage<T>>(std::vector<SurfDiscPtr>{surf}, name, func, comm);
}

template <typename T>
PostProcessorPtr makePostProcessorSurfaceIntegralAverage(const std::vector<SurfDiscPtr>& surfs, const std::string& name, T func, MPI_Comm comm)
{
  return std::make_shared<PostProcessorSurfaceIntegralAverage<T>>(surfs, name, func, comm);
}


class PostProcessorBCFlux : public PostProcessorBase
{
  public:

    PostProcessorBCFlux(const std::string& name, NeumannBCPtr bc, MPI_Comm comm) :
      m_name(name),
      m_bc(bc),
      m_comm(comm)
    {}

    int numValues() const override { return 1; }

    std::vector<std::string> getNames() const override { return {m_name}; }

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    std::string m_name;
    NeumannBCPtr m_bc;
    MPI_Comm m_comm;
};


class PostProcessorAirWindSkyBCFlux : public PostProcessorBCFlux
{
  public:

    PostProcessorAirWindSkyBCFlux(const std::string& name, std::shared_ptr<Heat::AirWindSkyNeumannBC> bc,
                                  Heat::HeatEquationSolar* heat_eqn_solar, MPI_Comm comm) :
      PostProcessorBCFlux(name, bc, comm),
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
                                  Heat::HeatEquationSolar* heat_eqn_solar, MPI_Comm comm) :
      m_name_prefix(name_prefix),
      m_bc(bc),
      m_heat_eqn_solar(heat_eqn_solar),
      m_comm(comm)
    {}

    // returns number of values this postprocessor returns
    int numValues() const override;

    std::vector<std::string> getNames() const override;

    std::vector<double> getValues(DiscVectorPtr u, AuxiliaryEquationsStoragePtr u_aux, double t) override;

  private:
    std::string m_name_prefix;
    std::shared_ptr<Heat::CombinedAirWindSkyNeumannBC> m_bc;
    Heat::HeatEquationSolar* m_heat_eqn_solar;    
    MPI_Comm m_comm;
};

}  // namespace

#endif