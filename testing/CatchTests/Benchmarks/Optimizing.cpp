#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <Benchmark.hpp>
#include <LaserTherm/HeatSolvers/_1D/Cartesian/CrankNicholson.hpp>

TEST_CASE("CrankNicholson Heat Solver Optimizations","[.][benchmarks]]")
{
  // this test case is for benchmarking the crank-nicholson
  // heat solver and comparing them
  BM::PerformanceBenchmark bm("reduce-copies");
  BM::Benchmark            meter;

  HeatSolvers::_1D::Cartesian::CrankNicholson<double> HeatSolver(1000);
  HeatSolver.T.setCoordinateSystem(Uniform(-5, 5));
  HeatSolver.VHC.set(4);
  HeatSolver.k.set(0.006);

  HeatSolver.T.set_f([&](auto i, auto cs) {
    auto x = cs->getCoord(i[0]);
    return 10 * exp(-2 * x * x);
  });

  double dt = HeatSolver.calcMaxTimeStep();

  meter.run([&]() { HeatSolver.stepForward(dt); });
  auto result = bm(meter);
  std::cout << "speedup over baseline: "
            << result.speedup_over_baseline.nominal() << "\n";
  std::cout << "speedup over  minimum: "
            << result.speedup_over_minimum.nominal() << "\n";
}
