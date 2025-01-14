#include "discretization/disc_vector.h"
#include "physics/heat/HeatEquationSolar.h"
#include "physics/heat/AuxiliaryEquationsSolar.h"
#include "physics/heat/basis_vals.h"
#include "physics/heat/source_terms_def.h"

namespace Heat {

HeatEquationSolar::HeatEquationSolar(DiscPtr disc, std::shared_ptr<SolarPositionCalculator> solar_position,
                  std::shared_ptr<EnvironmentInterface> environment_interface,
                  std::shared_ptr<InteriorAirTemperatureUpdator> air_temp_updator,
                  MPI_Comm comm)
: HeatEquation(disc, comm),
  m_solar_position(solar_position),
  m_environment(environment_interface),
  m_air_temp(air_temp_updator),
  m_aux_equations(std::make_shared<AuxiliaryEquationsSolar>(*this, air_temp_updator))
{}

//TODO: arguments are unused?
void HeatEquationSolar::initialize() 
{ 
  HeatEquation::initialize();
  m_air_temp->initialize(this); 
}


void HeatEquationSolar::addNeumannBC(NeumannBCPtr bc, bool is_exterior)
{
  addNeumannBC(bc); 
  m_is_neumann_bc_exterior.push_back(is_exterior);
}


void HeatEquationSolar::setTimeParameters(Real t, Real interior_air_temp)
{
  DirectionCosines solar_dir = m_solar_position->computePosition(t);

  EnvironmentData env_data   = m_environment->getEnvironmentData(t);
  const auto& neumann_bcs    = getNeumannBCs();
  for (size_t i=0; i < neumann_bcs.size(); ++i)
  {
    auto bc_air_wind_sky = std::dynamic_pointer_cast<AirWindSkyNeumannBC>(neumann_bcs[i]);
    if (bc_air_wind_sky)
    {
      if (m_is_neumann_bc_exterior[i])
      {
        bc_air_wind_sky->setAirTemperature(env_data.air_temp);
        bc_air_wind_sky->setAirSpeed(env_data.air_speed);
        bc_air_wind_sky->setAirDirection(env_data.air_direction);
        bc_air_wind_sky->setIRHorizontalRadiation(env_data.ir_horizontal_radiation);
        bc_air_wind_sky->setDirectNormalRadiation(env_data.direct_normal_radiation);
        bc_air_wind_sky->setDiffuseRadiation(env_data.diffuse_radiation);
        bc_air_wind_sky->setSolarDirection(solar_dir);
      } else
      {
        bc_air_wind_sky->setAirTemperature(interior_air_temp);
        bc_air_wind_sky->setAirSpeed(0);
        bc_air_wind_sky->setAirDirection(std::array<Real, 3>{1, 0, 0});
        bc_air_wind_sky->setIRHorizontalRadiation(0);
        bc_air_wind_sky->setDirectNormalRadiation(env_data.direct_normal_radiation);
        bc_air_wind_sky->setDiffuseRadiation(env_data.diffuse_radiation);
        bc_air_wind_sky->setSolarDirection(solar_dir);         
      }
    }
  }

  for (int i=0; i < getDiscretization()->getNumVolDiscs(); ++i)
  {
    if (hasSourceTerm(i))
    {
      SourceTermPtr src_term = getSourceTerm(i);
      auto src_term_air_wind_sky = std::dynamic_pointer_cast<SourceTermAirWindSky>(src_term);
      if (src_term_air_wind_sky)
      {
        src_term_air_wind_sky->setInteriorAirTemperature(interior_air_temp);
        src_term_air_wind_sky->setExteriorAirTemperature(env_data.air_temp);
        src_term_air_wind_sky->setAirSpeed(env_data.air_speed);
        src_term_air_wind_sky->setAirDirection(env_data.air_direction);
        src_term_air_wind_sky->setIRHorizontalRadiation(env_data.ir_horizontal_radiation);
        src_term_air_wind_sky->setDirectNormalRadiation(env_data.direct_normal_radiation);
        src_term_air_wind_sky->setDiffuseRadiation(env_data.diffuse_radiation);
        src_term_air_wind_sky->setSolarDirection(solar_dir);        
      }
    }
  }

  m_env_data = env_data;
  m_air_temp->setExteriorTemperature(env_data.air_temp);
}


AuxiliaryEquationsPtr HeatEquationSolar::getAuxEquations()
{ 
  return m_aux_equations; 
}

AuxiliaryEquationsSolarPtr HeatEquationSolar::getAuxEquationsSolar()
{ 
  return m_aux_equations; 
}


void HeatEquationSolar::computedRdTinterior_airProduct(DiscVectorPtr u, Real interior_temp, Real t, Real x, ArrayType<Real, 1>& b)
{
  //std::cout << "x = " << x << std::endl;
  //TODO: cache these
  setTimeParameters(t, interior_temp);
  auto rhs = makeDiscVector(getDiscretization());
  auto rhs_dot = makeDiscVector(getDiscretization());
  rhs_dot->set(0);

  computeNeumannBC_dotTair(*this, t, u, interior_temp, x, rhs, rhs_dot);
  computeSourceTerm_dotTair(*this, t, x, rhs, rhs_dot);
  rhs_dot->syncArrayToVector();

  auto& rhs_dot_vec = rhs_dot->getVector();
  assertAlways(b.shape()[0] == rhs_dot_vec.shape()[0], "vector sizes are incompatible");
  for (size_t i=0; i < rhs_dot_vec.shape()[0]; ++i)
  {
    //std::cout << "dof " << i << ", result = " << rhs_dot_vec[i] << std::endl;
    b[i] = rhs_dot_vec[i];
  }
}



void computeNeumannBC_dotTair(const HeatEquationSolar& physics, const Real t, DiscVectorPtr u, Real t_interior, Real t_interior_dot, 
                              DiscVectorPtr rhs, DiscVectorPtr rhs_dot)
{
  const auto& neumann_bcs = physics.getNeumannBCs();
  for (size_t i=0; i < neumann_bcs.size(); ++i)
    if (!physics.isNeumannBCExterior(i))
    {
      auto bc = neumann_bcs[i];
      computeNeumannBC_dotTair(bc, u, t_interior, t_interior_dot, t, rhs, rhs_dot);
    }

  rhs->markArrayModified();
  rhs_dot->markArrayModified();
}

void computeNeumannBC_dotTair(NeumannBCPtr bc, DiscVectorPtr u, Real t_interior, Real t_interior_dot, const Real t, 
                              DiscVectorPtr rhs, DiscVectorPtr rhs_dot)
{
  // Note: we currently don't require that a given surface has a single volume as
  //       as its upward adjacency, so we have to look up a different volume disc
  //       for every face
  
  auto bc_air_wind_sky = std::dynamic_pointer_cast<AirWindSkyNeumannBC>(bc);
  if (!bc_air_wind_sky)
    return;

  bc_air_wind_sky->setAirTemperature(t_interior);


  auto surf = bc->getSurfDisc();
  ArrayType<Real, 1> u_quad(boost::extents[surf->getNumQuadPtsPerFace()]);
  std::vector<Real> flux_vals(surf->getNumQuadPtsPerFace() * 3);
  std::vector<Real> flux_vals_dot(3*surf->getNumQuadPtsPerFace());

  Quadrature& quad = surf->quad;
  BasisVals2D basis(surf->face_group.getTPMapperSol(), quad.getPoints(), surf->face_group.getFaceNodesSol(), surf->face_group.ref_el_sol);

  for (int face=0; face < surf->getNumFaces(); ++face)
  {
    auto& face_spec   = surf->face_group.faces[face];
    auto& u_arr       = u->getArray(face_spec.vol_group);
    auto& res_arr     = rhs->getArray(face_spec.vol_group);
    auto& res_arr_dot = rhs_dot->getArray(face_spec.vol_group);
    std::fill(flux_vals_dot.begin(), flux_vals_dot.end(), 0); //TODO: should be unnecessary


    auto u_el = u_arr[boost::indices[face_spec.el_group][range()]];
    //auto res_el = res_arr[boost::indices[face_spec.el_group][range()]];

    surf->interp_vsq_flat[face_spec.face].interpolateVals(u_el, u_quad);
    bc_air_wind_sky->getValuedTair(face, t, u_quad.data(), flux_vals.data(), flux_vals_dot.data());


    for (int ki=0; ki < quad.getNumPoints(); ++ki)
      for (int kj=0; kj < quad.getNumPoints(); ++kj)
      {
        int k    = surf->quad_tp_nodemap[ki][kj];
        Real weight = quad.getWeight(ki) * quad.getWeight(kj);
        Real flux_normal = 0, flux_normal_dot = 0;
        for (int d=0; d < 3; ++d)
        {
          flux_normal     += surf->normals[face][k][d] * flux_vals[k + surf->getNumQuadPtsPerFace() * d];
          flux_normal_dot += surf->normals[face][k][d] * flux_vals_dot[k + surf->getNumQuadPtsPerFace() * d];
        }

        Real val = weight * flux_normal;
        Real val_dot = weight * flux_normal_dot * t_interior_dot;
        for (int i=0; i < surf->getNumSolPtsPerFace(); ++i)
        {
          int node_sol = surf->face_group.getFaceNodesSol()[face_spec.face][i];
          res_arr[face_spec.el_group][node_sol]     += basis.getValue(face_spec.face, i, ki, kj) * val;
          res_arr_dot[face_spec.el_group][node_sol] += basis.getValue(face_spec.face, i, ki, kj) * val_dot;
        }          
      }
  }
}


void computeSourceTerm_dotTair(const HeatEquation& physics, Real t, Real t_interior_dot,
                               DiscVectorPtr rhs, DiscVectorPtr rhs_dot)
{
  auto disc = physics.getDiscretization();

  for (int i=0; i < disc->getNumVolDiscs(); ++i)
  {
    if (physics.hasSourceTerm(i))
    {
      auto vol_disc = disc->getVolDisc(i);
      auto src_term = physics.getSourceTerm(i);
      auto src_term_air_wind_sky = std::dynamic_pointer_cast<SourceTermAirWindSky>(src_term);


      if (src_term_air_wind_sky)
      {
        auto& rhs_arr     = rhs->getArray(i);
        auto& rhs_arr_dot = rhs_dot->getArray(i);

        if (vol_disc->getNumElems() == 0)
          continue;

        computeSourceTerm_dotTair(vol_disc, src_term_air_wind_sky, t, t_interior_dot, rhs_arr, rhs_arr_dot);
      }
    }
  }

  rhs->markArrayModified();
  rhs_dot->markArrayModified();  
}

void computeSourceTerm_dotTair(const VolDiscPtr vol_disc, SourceTermAirWindSkyPtr src, Real t, Real t_interior_dot,
                               ArrayType<Real, 2>& rhs_arr, ArrayType<Real, 2>& rhs_arr_dot)
{
  auto& detJ  = vol_disc->detJ;
  BasisVals basis_vals(vol_disc->vol_group.getTPMapperSol(), Mesh::TensorProductMapper(vol_disc->quad.getPoints()));
  std::vector<Real> src_vals(vol_disc->getNumQuadPtsPerElement()), src_vals_dot(vol_disc->getNumQuadPtsPerElement());

  auto& rev_nodemap = basis_vals.getRevNodemapOut();
  for (int el=0; el < vol_disc->getNumElems(); ++el)
  {
    src->getValues(el, t, src_vals.data());
    src->getValues_dTair(el, t, t_interior_dot, src_vals_dot.data());

    for (int q=0; q < vol_disc->getNumQuadPtsPerElement(); ++q)
    {
      int k_i = rev_nodemap[q][0]; int k_j = rev_nodemap[q][1]; int k_k = rev_nodemap[q][2];
      Real weight = vol_disc->quad.getWeight(k_i) * vol_disc->quad.getWeight(k_j) * vol_disc->quad.getWeight(k_k);

      for (int i=0; i < vol_disc->getNumSolPtsPerElement(); ++i)
      {
        Real basis_val      =  basis_vals.getValue(i, q);
        rhs_arr[el][i]     += basis_val * weight * src_vals[q] / detJ[el][q];
        rhs_arr_dot[el][i] += basis_val * weight * src_vals_dot[q] / detJ[el][q];
      }
    }
  }
}

}