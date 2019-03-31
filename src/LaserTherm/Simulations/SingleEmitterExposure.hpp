#ifndef LaserTherm_Simulations_Basic_hpp
#define LaserTherm_Simulations_Basic_hpp

/** @file Basic.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/02/19
 */

#include <iomanip>
#include <iostream>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/signals2/signal.hpp>

#include "../Utils/TypeTraits.hpp"

namespace Simulations
{

template<class EmitterType, class HeatSolverType>
class SingleEmitterExposure
{
 public:
  using real_type = typename get_real_type<decltype(HeatSolverType::A)>::type;

  EmitterType    emitter;
  HeatSolverType heat_solver;

  // time accumulation should be done in a minimum of double precision
  // even if thermal solver is in float percision
  // otherwise we will accumulate significant error
  // see https://www.youtube.com/watch?v=k12BJGSc2Nc&feature=youtu.be&t=1513
  long double tmax;
  long double dt;
  long double t;

  void run();

  boost::signals2::signal<void(const decltype(HeatSolverType::T)&, const long double&)> sig_broadcastTemperatureBeforeHeatSolverStepForward;
  boost::signals2::signal<void(const decltype(HeatSolverType::A)&, const long double&)> sig_broadcastSourceTermBeforeHeatSolverStepForward;
};

template<class EmitterType, class HeatSolverType>
void SingleEmitterExposure<EmitterType, HeatSolverType>::run()
{
  uint64_t Nt = 1 + tmax / dt;

  t = 0;
  for (int i = 0; i < Nt; ++i) {
    t += dt;
    sig_broadcastTemperatureBeforeHeatSolverStepForward( heat_solver.T, t );
    sig_broadcastSourceTermBeforeHeatSolverStepForward( heat_solver.A, t );
    heat_solver.stepForward(dt);
  }
}

}  // namespace Simulations

#endif  // include protector
