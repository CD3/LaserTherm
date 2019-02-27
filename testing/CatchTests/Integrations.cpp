#include "catch.hpp"

#include<iostream>
#include<fstream>

#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>
#include <LaserTherm/HeatSources/_1D/BeersLaw.hpp>


TEST_CASE("Crank-Nicholson Heat Solver with Beers Law Source")
{

  HeatSolvers::_1D::CrankNicholson<double> HeatSolver(600);
  HeatSolver.T.setCoordinateSystem( Uniform(-2,2) );

  HeatSources::_1D::BeersLaw<double> HeatSource(600);
  HeatSource.A.setCoordinateSystem( Uniform(-2,2) );


  HeatSource.E0 = 1;
  HeatSource.mu_a.set_f(
      [](auto i, auto cs)
      {
      if(i[0] < 300)
        return 0;
      return 100;
      }
      );
  HeatSolver.VHC.set( 4.1868 );
  HeatSolver.k.set( 0.00628 );

  HeatSource.compute();
  HeatSolver.A = HeatSource.A;


  std::ofstream out("T.txt");
  out << HeatSolver.T;
  double dt = HeatSolver.calcMaxTimeStep();
  for(int i = 0; i < 1000; ++i)
  {
    HeatSolver.stepForward(dt);
    out << "\n" << HeatSolver.T;
  }

}
