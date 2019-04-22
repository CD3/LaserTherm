#ifndef LaserTherm_Configuration_Builders_hpp
#define LaserTherm_Configuration_Builders_hpp

/** @file Builders.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/24/19
 */

#include <iostream>
#include <map>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "../Configuration/Manager.hpp"
#include "../Emitters/Basic.hpp"
#include "../HeatSolvers/_1D/CrankNicholson.hpp"
#include "../MaterialStructure.hpp"
#include "../Materials/Basic.hpp"
#include "../Simulations/SingleEmitterExposure.hpp"
#include "../Structures/_1D/AnyStructure.hpp"
#include "../Structures/_1D/Infinite.hpp"
#include "../Structures/_1D/Slab.hpp"
#include "../Waveforms/ContinuousWave.hpp"

namespace Builders
{
template<class REAL>
void build(Materials::Basic<REAL>&            material,
           const boost::property_tree::ptree& config)
{
  boost::optional<REAL> prop;

  /* see ../Materials/Basic.hpp for list of properties */
  prop = config.get_optional<REAL>("thermal.specific_heat");
  if (prop) material.setSpecificHeatCapacity(prop.value());

  prop = config.get_optional<REAL>("thermal.density");
  if (prop) material.setDensity(prop.value());

  prop = config.get_optional<REAL>("thermal.conductivity");
  if (prop) material.setThermalConductivity(prop.value());

  prop = config.get_optional<REAL>("optical.absorption_coefficient");
  if (prop) material.setAbsorptionCoefficient(prop.value());

  prop = config.get_optional<REAL>("optical.scattering_coefficient");
  if (prop) material.setScatteringCoefficient(prop.value());
}

template<class REAL>
void build(std::map<std::string, Materials::Basic<REAL>>& materials,
           const boost::property_tree::ptree&             config)
{
  using MaterialType = Materials::Basic<REAL>;

  for (auto it : config) {
    MaterialType mat;
    build(mat, it.second);
    materials[it.first] = mat;
  }
}

/**
 * Builds an vector of 1D material layers from a configuration defined in a
 * boost property_tree.
 */
template<class REAL>
void build(std::vector<MaterialStructure<Materials::Basic<REAL>,
                                         Structures::_1D::AnyStructure<REAL>>>&
                                              structures,
           const boost::property_tree::ptree& config)
{
  using MaterialType          = Materials::Basic<REAL>;
  using StructureType         = Structures::_1D::AnyStructure<REAL>;
  using MaterialStructureType = MaterialStructure<MaterialType, StructureType>;
  using InfStucture           = Structures::_1D::Infinite<REAL>;
  using LayerStucture         = Structures::_1D::Slab<REAL>;

  std::map<std::string, MaterialType> materials;
  auto materials_config = config.get_child_optional("materials");
  if (materials_config) {
    build(materials, materials_config.value());
  }

  auto layers_config = config.get_child_optional("layers");
  if (layers_config) {
    boost::optional<REAL> current_position;
    for (auto it : layers_config.value()) {
      MaterialType mat;

      if (it.second.find("material") != it.second.not_found())
        mat = materials[it.second.get<std::string>("material")];
      else
        build(mat, it.second);

      auto position  = it.second.get_optional<REAL>("position");
      auto thickness = it.second.get_optional<REAL>("thickness");

      if (!position && !thickness && !current_position) {
        structures.push_back(
            MaterialStructureType(mat, StructureType(InfStucture())));
        continue;
      }

      LayerStucture layer;
      if (!position && current_position) {
        layer.setMinSurfacePosition(current_position.value());
      }

      if (position && !current_position) {
        layer.setMinSurfacePosition(position.value());
      }

      if (thickness) {
        layer.setThickness(thickness.value());
        current_position = layer.getMaxSurfacePosition();
      }

      structures.push_back(MaterialStructureType(mat, layer));
    }
  }
}

template<class EmitterType, class HeatSolverType>
void build(Simulations::SingleEmitterExposure<EmitterType, HeatSolverType>& sim,
           const Configuration::Manager& config)
{
  using real_type = typename get_real_type<HeatSolverType>::type;

  sim.tmax = config.get<long double>("simulation.time.end");
  sim.dt   = config.get<long double>("simulation.time.dt.max");

  sim.heat_solver =
      decltype(sim.heat_solver)(config.get<size_t>("simulation.grid.x.n"));
  if (config.get<std::string>("simulation.grid.type") == "uniform") {
    sim.heat_solver.T.setCoordinateSystem(
        Uniform(config.get_quantity<real_type>("simulation.grid.x.min", "cm"),
                config.get_quantity<real_type>("simulation.grid.x.max", "cm")));
  }

  sim.emitter.A.reset(sim.heat_solver.T.getCoordinateSystemPtr());
  sim.emitter.mu_a.reset(sim.heat_solver.T.getCoordinateSystemPtr());

  sim.heat_solver.A.set(0);
  sim.emitter.A.set(0);

  // sim.emitter.setIrradiance(config.get<real_type>("emitters.0.irradiance"));
  // sim.emitter.setStartTime(config.get<real_type>("emitters.0.start"));
  // sim.emitter.setExposureDuration(
  // config.get<real_type>("emitters.0.exposure_duration"));
}

}  // namespace Builders

#endif  // include protector
