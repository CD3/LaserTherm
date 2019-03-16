#ifndef LaserTherm_Simulations_Basic_hpp
#define LaserTherm_Simulations_Basic_hpp

/** @file Basic.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/02/19
 */

#include <iostream>
#include <iomanip>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "../Utils/TypeTraits.hpp"

namespace Simulations
{
using conf_tree = boost::property_tree::ptree;

template<class EmitterType, class HeatSolverType>
class SingleEmitterExposure
{
 public:
  using real_type = typename get_real_type<decltype(HeatSolverType::A)>::type;

  conf_tree      config;
  EmitterType    emitter;
  HeatSolverType heat_solver;

  void configure(const conf_tree& config);
  void run();
};

template<class EmitterType, class HeatSolverType>
void SingleEmitterExposure<EmitterType, HeatSolverType>::configure(
    const conf_tree& a_config)
{
  this->config = a_config;

  int n = config.get<int>("simulation.grid.x.n");

  heat_solver.reset(n);

  real_type xmin = config.get<real_type>("simulation.grid.x.min");
  real_type xmax = config.get<real_type>("simulation.grid.x.max");

  heat_solver.T.setCoordinateSystem(Uniform(xmin, xmax));

  emitter.A.reset(heat_solver.T.getCoordinateSystemPtr());
  emitter.mu_a.reset(heat_solver.T.getCoordinateSystemPtr());
  heat_solver.A.reset(heat_solver.T.getCoordinateSystemPtr());
  heat_solver.VHC.reset(heat_solver.T.getCoordinateSystemPtr());
  heat_solver.k.reset(heat_solver.T.getCoordinateSystemPtr());

  // configure layers
  //structure.set<SpecificHeat>( heat_solver.k );
  //structure.set<VolumetricHeatCapacity>( heat_solver.VHC );
  //structure.set<AbsorptionCoefficient>( emitter.mu_a );

  // configure boundary conditions...

  // connect emitters computeSourceTerm function to the heat solvers
  // sig_askSourceTerm signal

  heat_solver.sig_askSourceTerm.connect(
      [](auto& A, const auto& t) { A.set(0); });
  heat_solver.sig_askSourceTerm.connect(
      [this](auto& A, const auto& t) { this->emitter.addSourceTerm(A, t); });
}

template<class EmitterType, class HeatSolverType>
void SingleEmitterExposure<EmitterType, HeatSolverType>::run()
{
  // time accumulation should be done in a minimum of double precision
  // even if thermal solver is in float percision
  // otherwise we will accumulate significant error
  // see https://www.youtube.com/watch?v=k12BJGSc2Nc&feature=youtu.be&t=1513
  long double tmax = config.get<long double>("simulation.time.end");
  long double dt = config.get<long double>("simulation.time.dt.min");
  uint64_t Nt = 1 + tmax / dt;

  double t = 0;
  for(int i = 0; i < Nt; ++i)
  {
    t += dt;

    heat_solver.stepForward(dt);
  }

}

}  // namespace Simulations

#endif  // include protector
