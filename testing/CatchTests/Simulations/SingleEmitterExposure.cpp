#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <LaserTherm/Configuration/Builders.hpp>
#include <LaserTherm/Configuration/Manager.hpp>
#include <LaserTherm/Emitters/Basic.hpp>
#include <LaserTherm/HeatSolvers/_1D/Cartesian/CrankNicholson.hpp>
#include <LaserTherm/HeatSources/_1D/Cartesian/BeersLaw.hpp>
#include <LaserTherm/MaterialStructure.hpp>
#include <LaserTherm/Materials/Basic.hpp>
#include <LaserTherm/Simulations/SingleEmitterExposure.hpp>
#include <LaserTherm/Structures/_1D/AnyStructure.hpp>
#include <LaserTherm/Structures/_1D/Cartesian/Slab.hpp>
#include <LaserTherm/Structures/_1D/Infinite.hpp>
#include <LaserTherm/Waveforms/ContinuousWave.hpp>

TEST_CASE("Simple Simulation Test", "[.][simulations][longrunning]") {
  // this is basically a testing ground for building and running a simulation

  Configuration::Manager config;

  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Length>("cm");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Mass>("g");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Time>("s");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Temperature>(
      "K");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Amount>("mol");
  config.unit_registry
      .addBaseUnit<UnitConvert::Dimension::Name::ElectricalCurrent>("A");
  config.unit_registry
      .addBaseUnit<UnitConvert::Dimension::Name::LuminousIntensity>("cd");

  config.unit_registry.addUnit("m = 100 cm");
  config.unit_registry.addUnit("L = 1000 cm^3");
  config.unit_registry.addUnit("J = kg m^2 / s^2");
  config.unit_registry.addUnit("W = J/s");
  config.unit_registry.addUnit("cal = 4.184 J");
  config.unit_registry.addUnit("degC = K - 273.15");
  config.unit_registry.addUnit("degK = K");

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
  materials.water.thermal.conductivity = 0.00628 W/cm/K
  materials.water.thermal.specific_heat = 4.1868 J/g/degK
  materials.water.thermal.bc.convection.transfer_rate = 1e-3 W/cm^s/K

  layers.1.description = absorber
  layers.1.optical.absorption_coefficient = 1 cm^-1
  layers.1.thermal.conductivity = 1 J/g/degK
  layers.1.thermal.density = 2 g/mL
  layers.1.thermal.specific_heat = 3 J/g/degK
  layers.1.position = 0 cm
  layers.1.thickness = 5 mm
  )";
  std::stringstream in(config_text);

  std::ofstream out("config.ini");
  out << config_text;
  out.close();

  config.load("config.ini");
  convertPropertyTreeUnits(config.configuration, config.unit_registry);

  Simulations::SingleEmitterExposure<
      Emitters::Basic<HeatSources::_1D::Cartesian::BeersLaw<double>,
                      Waveforms::ContinuousWave<double>>,
      HeatSolvers::_1D::Cartesian::CrankNicholson<double>>
      sim;

  size_t xN = config.get<size_t>("simulation.grid.x.n");
  double xmin = config.get<double>("simulation.grid.x.min");
  double xmax = config.get<double>("simulation.grid.x.max");
  sim.tmax = config.get<long double>("simulation.time.end");
  sim.dt = config.get<long double>("simulation.time.dt.max");

  sim.heat_solver =
      decltype(sim.heat_solver)(config.get<size_t>("simulation.grid.x.n"));
  sim.heat_solver.T.setCoordinateSystem(
      Uniform(config.get<double>("simulation.grid.x.min"),
              config.get<double>("simulation.grid.x.max")));

  sim.emitter.A.reset(sim.heat_solver.T.getCoordinateSystemPtr());
  sim.emitter.mu_a.reset(sim.heat_solver.T.getCoordinateSystemPtr());

  sim.heat_solver.A.set(0);
  sim.emitter.A.set(0);

  sim.emitter.setIrradiance(config.get<double>("emitters.0.irradiance"));
  sim.emitter.setStartTime(config.get<double>("emitters.0.start"));
  sim.emitter.setExposureDuration(
      config.get<double>("emitters.0.exposure_duration"));

  using MaterialType = Materials::Basic<double>;
  std::vector<
      MaterialStructure<MaterialType, Structures::_1D::AnyStructure<double>>>
      structures;

  Builders::build(structures, config.configuration);
  CHECK(structures.size() == 2);

  for (auto &s : structures) {
    // set conductivity
    sim.heat_solver.k.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getThermalConductivity())
        return s.material.getThermalConductivity().value();
      return std::nullopt;
    });

    // set volumetric heat capacity
    sim.heat_solver.VHC.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0])) {
        if (s.material.getDensity() && s.material.getSpecificHeatCapacity()) {
          return s.material.getDensity().value() *
                 s.material.getSpecificHeatCapacity().value();
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
  CHECK(sim.heat_solver.k(0) == Approx(1 * 1000 * 100 * 100));
  CHECK(sim.heat_solver.VHC(0) == Approx(6 * 1000 * 100 * 100));
  CHECK(sim.emitter.mu_a(0) == Approx(1));

  sim.heat_solver.T.set_f([](auto ind, auto cs) {
    auto x = cs->getCoord(ind[0]);
    return 10 * exp(-x * x / 5);
  });

  sim.heat_solver.sig_askSourceTerm.connect(
      [](auto &A, const auto &t) { A.set(0); });

  sim.heat_solver.sig_askSourceTerm.connect([&sim](auto &A, const auto &t) {
    sim.emitter.addSourceTerm(A, sim.t + t);
  });

  std::ofstream T_out("Simulation-T.txt");
  std::ofstream A_out("Simulation-A.txt");

  sim.sig_broadcastTemperatureBeforeHeatSolverStepForward.connect(
      [&T_out](auto T, auto t) {
        static long double last_time = -1;
        if (last_time < 0 || t - last_time > 0.5) {
          T_out << "# t = " << t << "\n";
          T_out << T;
          T_out << "\n";
          last_time = t;
        }
      });
  sim.sig_broadcastSourceTermBeforeHeatSolverStepForward.connect(
      [&A_out](auto A, auto t) {
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

TEST_CASE("Simulation Builder Test") {
  // this is basically a testing ground for building and running a simulation

  Configuration::Manager config;

  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Length>("cm");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Mass>("g");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Time>("s");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Temperature>(
      "K");
  config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Amount>("mol");
  config.unit_registry
      .addBaseUnit<UnitConvert::Dimension::Name::ElectricalCurrent>("A");
  config.unit_registry
      .addBaseUnit<UnitConvert::Dimension::Name::LuminousIntensity>("cd");

  config.unit_registry.addUnit("m = 100 cm");
  config.unit_registry.addUnit("L = 1000 cm^3");
  config.unit_registry.addUnit("J = kg m^2 / s^2");
  config.unit_registry.addUnit("W = J/s");
  config.unit_registry.addUnit("cal = 4.184 J");
  config.unit_registry.addUnit("degC = K - 273.15");

  std::string config_text = R"(
  simulation.dimensions = 1

  simulation.grid.type = "uniform"
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
  materials.water.thermal.conductivity = 0.00628 W/cm/K
  materials.water.thermal.specific_heat = 4.1868 J/g/degK
  materials.water.thermal.bc.convection.transfer_rate = 1e-3 W/cm^s/K

  layers.1.description = absorber
  layers.1.optical.absorption_coefficient = 1 cm^-1
  layers.1.thermal.conductivity = 1 W/cm/K
  layers.1.thermal.density = 2 g/mL
  layers.1.thermal.specific_heat = 3 J/g/degK
  layers.1.position = 0 mm
  layers.1.thickness = 5 mm
  )";
  std::stringstream in(config_text);

  std::ofstream out("config.ini");
  out << config_text;
  out.close();

  config.load("config.ini");
  convertPropertyTreeUnits(config.configuration, config.unit_registry);

  Simulations::SingleEmitterExposure<
      Emitters::Basic<HeatSources::_1D::Cartesian::BeersLaw<double>,
                      Waveforms::ContinuousWave<double>>,
      HeatSolvers::_1D::Cartesian::CrankNicholson<double>>
      sim;

  Builders::build(sim, config);
}
