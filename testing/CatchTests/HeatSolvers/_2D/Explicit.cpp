#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <LaserTherm/HeatSolvers/_2D/Cylindrical/Explicit.hpp>

TEST_CASE("Explicit 2D Cylindrical Heat Solver")
{

     //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
    Explicit<double> HeatSolver(10,20);
/*     HeatSolver.T.setCoordinateSystem(Uniform(0, 2),Uniform(0,4)); */

/*     HeatSolver.T.set(1.0); */
/*     HeatSolver.A.set(2.0); */
/*     HeatSolver.VHC.set(3.0); */
/*     HeatSolver.k.set(4.0); */


/*     HeatSolver.stepForward(1); */

}
