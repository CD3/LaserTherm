#include "catch.hpp"

#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>
#include <LaserTherm/HeatSources/_1D/BeersLaw.hpp>


TEST_CASE("Crank-Nicholson Heat Solver with Beers Law Source")
{

  HeatSolvers::_1D::CrankNicholson<double> HeatSolver(11);
  HeatSolver.T.setCoordinateSystem( Uniform(-2,2) );

  HeatSources::_1D::BeersLaw<double> HeatSource(11);
  HeatSource.A.setCoordinateSystem( Uniform(-2,2) );


}
