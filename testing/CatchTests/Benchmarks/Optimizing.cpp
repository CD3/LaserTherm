#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>

#include <Benchmark.hpp>
#include <LaserTherm/HeatSolvers/_1D/Cartesian/CrankNicholson.hpp>
#include <LaserTherm/HeatSolvers/_2D/Cylindrical/UniformExplicit.hpp>

TEST_CASE("CrankNicholson Heat Solver Optimizations", "[.][benchmarks]") {
  // this test case is for benchmarking the crank-nicholson
  // heat solver and comparing them
  BM::PerformanceBenchmark bm("reduce-copies");
  BM::Benchmark meter;

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

TEST_CASE("2D Cylindrical Explicit Heat Solver Optimizations",
          "[.][benchmarks]") {
  // this test case is for benchmarking the crank-nicholson
  // heat solver and comparing them
  BM::PerformanceBenchmark bm("2d-cyl-exp-stepforward");
  BM::Benchmark meter;

  double R = 2;
  double L = 4;
  int rN = 141;
  int zN = 282;
  double dR = R / rN;
  double dL = L / zN;
  UniformExplicit<double> HeatSolver(zN - 2, rN - 2);
  Field<double, 2> Aplot(zN, rN);
  HeatSolver.T.setCoordinateSystem(Uniform(dL, L - dL), Uniform(dR, R - dR));
  Aplot.setCoordinateSystem(Uniform(0., L), Uniform(0., R));

  double k = 0.5;    // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18;   // J/kg/K
  HeatSolver.A.set(0.0);
  HeatSolver.VHC.set(rho * c);
  HeatSolver.k.set(k);
  HeatSolver.minZBC.f = 0;
  HeatSolver.maxZBC.f = 0;
  HeatSolver.maxRBC.f = 0;

  BENCHMARK("stepForward") {
    HeatSolver.stepForward(0.01);
    return HeatSolver.T(0, 0);
  };

  meter.run([&]() { HeatSolver.stepForward(0.01); });
  auto result = bm(meter);
  std::cout << "speedup over baseline: "
            << result.speedup_over_baseline.nominal() << "\n";
  std::cout << "speedup over  minimum: "
            << result.speedup_over_minimum.nominal() << "\n";
}
