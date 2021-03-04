#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <LaserTherm/HeatSolvers/_2D/Cylindrical/Explicit.hpp>
using namespace std;

TEST_CASE("Explicit 2D Cylindrical Heat Solver")
{

     //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
    Explicit<double> HeatSolver(10,20);
    HeatSolver.T.setCoordinateSystem(Uniform(0, 2),Uniform(0,4));
    HeatSolver.A.setCoordinateSystem(Uniform(0, 2),Uniform(0,4));
    HeatSolver.VHC.setCoordinateSystem(Uniform(0, 2),Uniform(0,4));
    HeatSolver.k.setCoordinateSystem(Uniform(0, 2),Uniform(0,4));
//      Herbie.destination.address.directions("tomanek");

    HeatSolver.T.set(1.0);
    HeatSolver.A.set(2.0);
    HeatSolver.VHC.set(3.0);
    HeatSolver.k.set(4.0);

    {
      ofstream output;
      output.open("Tbefore.txt");
      output << HeatSolver.T;
      output.close();
    }
    
    HeatSolver.stepForward(1);
    
    {
      ofstream output;
      output.open("Tafter.txt");
      output << HeatSolver.T;
      output.close();
    }
}
