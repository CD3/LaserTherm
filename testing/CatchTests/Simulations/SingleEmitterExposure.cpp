#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/info_parser.hpp>

#include <LaserTherm/Emitters/Basic.hpp>
#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>
#include <LaserTherm/HeatSources/_1D/BeersLaw.hpp>
#include <LaserTherm/MaterialStructure.hpp>
#include <LaserTherm/Materials/Basic.hpp>
#include <LaserTherm/Simulations/SingleEmitterExposure.hpp>
#include <LaserTherm/Structures/_1D/AnyStructure.hpp>
#include <LaserTherm/Structures/_1D/Infinite.hpp>
#include <LaserTherm/Structures/_1D/Slab.hpp>
#include <LaserTherm/Waveforms/ContinuousWave.hpp>
#include <LaserTherm/Configuration/Builders.hpp>
#include <LaserTherm/Configuration/ptree_utils.hpp>

TEST_CASE("Simple Simulation Test")
{
  // this is basically a testing ground for building and running a simulation
  UnitRegistry ureg;
  ureg.addBaseUnit<Dimension::Name::Length>("cm");
  ureg.addBaseUnit<Dimension::Name::Mass>("g");
  ureg.addBaseUnit<Dimension::Name::Time>("s");
  ureg.addBaseUnit<Dimension::Name::Temperature>("K");
  ureg.addBaseUnit<Dimension::Name::Amount>("mol");
  ureg.addBaseUnit<Dimension::Name::ElectricalCurrent>("A");
  ureg.addBaseUnit<Dimension::Name::LuminousIntensity>("cd");

  ureg.addUnit("m = 100 cm");
  ureg.addUnit("L = 1000 cm^3");
  ureg.addUnit("J = kg m^2 / s^2");
  ureg.addUnit("W = J/s");
  ureg.addUnit("cal = 4.184 J");
  ureg.addUnit("degC = K - 273.15");

  boost::property_tree::ptree config;

  std::string config_text = R"(
  simulation.dimensions = 1

  simulation.grid.x.n = 100
  simulation.grid.x.min = 0 cm
  simulation.grid.x.max = 2 cm
  simulation.time.end = 3 s
  simulation.time.dt.max = 1e-3 s
  simulation.baseline_temperature = 37 degC

  simulation.boundary_condition.min.type = surface
  simulation.boundary_condition.max.type = sink

  emitters.0.irradiance = 5 W/cm^2
  emitters.0.wavelength = 1064 nm
  emitters.0.start = 0 s
  emitters.0.exposure_duration = 1 s

  layers.0.material = water

  materials.water.thermal.density = 1 g/mL
  materials.water.thermal.conductivity = 0.00628 W/cm^2/degK
  materials.water.thermal.specific_heat = 4.1868 J/g/degK
  materials.water.thermal.bc.convection.transfer_rate = 1e-3 W/cm^s/K

  layers.1.description = absorber
  layers.1.thermal.absorption_coefficient = 1 g/mL
  layers.1.position = 0 cm
  layers.1.thickness = 10e-4 cm
  )";
  std::stringstream in(config_text);

  boost::property_tree::read_ini( in, config );
  config = unflatten_ptree(config, '.', '|');
  convertPropertyTreeUnits( config, ureg );
  boost::property_tree::write_info( std::cout, config );
  



  Simulations::SingleEmitterExposure<
      Emitters::Basic<HeatSources::_1D::BeersLaw<double>,
                      Waveforms::ContinuousWave<double> >,
      HeatSolvers::_1D::CrankNicholson<double> >
      sim;

  size_t xN   = config.get<size_t>("simulation.grid.x.n");
  double xmin = config.get<double>("simulation.grid.x.min");
  double xmax = config.get<double>("simulation.grid.x.max");
  sim.tmax    = config.get<long double>("simulation.time.end");
  sim.dt      = config.get<long double>("simulation.time.dt.max");

  sim.heat_solver =
      decltype(sim.heat_solver)(config.get<size_t>("simulation.grid.x.n"));
  sim.heat_solver.T.setCoordinateSystem(
      Uniform(config.get<double>("simulation.grid.x.min"),
              config.get<double>("simulation.grid.x.max")));

  sim.emitter.A.reset(sim.heat_solver.T.getCoordinateSystemPtr());
  sim.emitter.mu_a.reset(sim.heat_solver.T.getCoordinateSystemPtr());

  sim.heat_solver.A.set(0);
  sim.emitter.A.set(0);

  using MaterialType = Materials::Basic<double>;
  std::vector<
      MaterialStructure<MaterialType, Structures::_1D::AnyStructure<double> > >
      structures;

  Builders::buildStructures(structures, config);

  for (auto& s : structures) {
    // set conductivity
    sim.heat_solver.k.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getThermalConductivity())
        return s.material.getThermalConductivity().value();
      return std::nullopt;
    });

    // set volumetric heat capacity
    sim.heat_solver.VHC.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0])) {
      if(s.material.getDensity() && s.material.getSpecificHeatCapacity())
      {
      return s.material.getDensity().value()*s.material.getSpecificHeatCapacity().value();
      }
      }

      return std::nullopt;
    });

    // set absorption coefficient
    sim.emitter.mu_a.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getAbsorptionCoefficient())
        return s.material.getAbsorptionCoefficient().value();
      return std::nullopt;
    });
  }

  sim.heat_solver.T.set_f([](auto ind, auto cs) {
    auto x = cs->getCoord(ind[0]);
    return 10 * exp(-x * x / 5);
  });

  std::ofstream T_out("Simulation-T.txt");
  std::ofstream A_out("Simulation-A.txt");

  sim.sig_broadcastTemperatureBeforeStep.connect([&T_out](auto T, auto t) {
    static long double last_time = -1;
    if (last_time < 0 || t - last_time > 0.5) {
      T_out << "# t = " << t << "\n";
      T_out << T;
      T_out << "\n";
      last_time = t;
    }
  });
  sim.sig_broadcastTemperatureBeforeStep.connect([&A_out](auto A, auto t) {
    static long double last_time = -1;
    if (last_time < 0 || t - last_time > 1) {
      A_out << "# t = " << t << "\n";
      A_out << A;
      A_out << "\n";
      last_time = t;
    }
  });

  sim.run();
}
