#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <boost/property_tree/ptree.hpp>

#include <LaserTherm/Emitters/Basic.hpp>
#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>
#include <LaserTherm/HeatSources/_1D/BeersLaw.hpp>
#include <LaserTherm/Simulations/SingleEmitterExposure.hpp>
#include <LaserTherm/Waveforms/ContinuousWave.hpp>

TEST_CASE("Simple Simulation Test")
{
  Simulations::SingleEmitterExposure<
      Emitters::Basic<HeatSources::_1D::BeersLaw<double>,
                      Waveforms::ContinuousWave<double> >,
      HeatSolvers::_1D::CrankNicholson<double> >
      sim;

  boost::property_tree::ptree config;

  config.put("simulation.dimensions", 1);

  config.put("simulation.grid.x.n", 100);
  config.put("simulation.grid.x.min", 0 /*cm*/);
  config.put("simulation.grid.x.max", 2 /*cm*/);
  config.put("simulation.time.end", 3 /*s*/);
  config.put("simulation.time.dt.max", 1e-3 /*s*/);
  config.put("simulation.baseline_temperature", 37 /*degC*/);

  config.put("simulation.boundary_condition.min.type", "surface");
  config.put("simulation.boundary_condition.max.type", "sink");

  config.put("emitters.0.irradiance", 5 /*W/cm^2*/);
  config.put("emitters.0.wavelength", 1064 /*nm*/);
  config.put("emitters.0.start", 0 /*s*/);
  config.put("emitters.0.exposure_duration", 1 /*s*/);

  config.put("layers.0.optical.absorption_coefficient", 50 /*1/cm*/);
  config.put("layers.0.thermal.density", 1 /*g/mL*/);
  config.put("layers.0.thermal.conductivity", 0.00628 /*W/cm^2/degK*/);
  config.put("layers.0.thermal.specific_heat", 4.1868 /*J/g/degK*/);
  config.put("layers.0.thermal.bc.convection.transfer_rate", 1e-3 /*W/cm^s/degK*/);


  size_t xN = config.get<size_t>("simulation.grid.x.n");
  double xmin = config.get<double>("simulation.grid.x.min");
  double xmax = config.get<double>("simulation.grid.x.max");
  sim.tmax = config.get<long double>("simulation.time.end");
  sim.dt = config.get<long double>("simulation.time.dt.max");

  sim.heat_solver = decltype(sim.heat_solver)( config.get<size_t>("simulation.grid.x.n")  );
  sim.heat_solver.T.setCoordinateSystem(Uniform(config.get<double>("simulation.grid.x.min"), config.get<double>("simulation.grid.x.max")));



  sim.emitter.A.reset(sim.heat_solver.T.getCoordinateSystemPtr());
  sim.emitter.mu_a.reset(sim.heat_solver.T.getCoordinateSystemPtr());


  sim.heat_solver.A.set(0);
  sim.heat_solver.k.set(1);
  sim.heat_solver.VHC.set(1);
  sim.emitter.A.set(0);

  sim.heat_solver.T.set_f( [](auto ind, auto cs)
      {
      auto x = cs->getCoord(ind[0]);
      return 10*exp(-x*x/5);
      });

  std::ofstream T_out("Simulation-T.txt");
  std::ofstream A_out("Simulation-A.txt");

  sim.sig_broadcastTemperatureBeforeStep.connect([&T_out](auto T, auto t)
      {
      static long double last_time = -1;
      if( last_time < 0 || t - last_time > 0.5)
      {
        T_out << "# t = " << t << "\n";
        T_out << T;
        T_out << "\n";
        last_time = t;
      }
      }
      );
  sim.sig_broadcastTemperatureBeforeStep.connect([&A_out](auto A, auto t)
      {
      static long double last_time = -1;
      if( last_time < 0 || t - last_time > 1)
      {
        A_out << "# t = " << t << "\n";
        A_out << A;
        A_out << "\n";
        last_time = t;
      }
      }
      );


  sim.run();
}
