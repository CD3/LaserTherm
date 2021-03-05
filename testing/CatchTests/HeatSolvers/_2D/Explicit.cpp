#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <LaserTherm/HeatSolvers/_2D/Cylindrical/Explicit.hpp>
using namespace std;

TEST_CASE("Explicit 2D Cylindrical Heat Solver")
{

     //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
    Explicit<double> HeatSolver(10,20);
    Field<double,2> Tinit(10,20);
    Tinit.setCoordinateSystem(Uniform(-2, 2),Uniform(-4,4));
    HeatSolver.T.setCoordinateSystem(Uniform(-2, 2),Uniform(-4,4));
    HeatSolver.A.setCoordinateSystem(Uniform(-2, 2),Uniform(-4,4));
    HeatSolver.VHC.setCoordinateSystem(Uniform(-2, 2),Uniform(-4,4));
    HeatSolver.k.setCoordinateSystem(Uniform(-2, 2),Uniform(-4,4));
//      Herbie.destination.address.directions("tomanek");
    auto Tlin = [](double z, double r) -> double {
      return exp(-(pow(r, 2) + pow(z, 2)));
    };

    HeatSolver.T.set_f([&Tlin] (auto v) -> double {return Tlin(v[0], v[1]);});
    Tinit.set_f([&Tlin] (auto v) -> double {return Tlin(v[0], v[1]);});
    HeatSolver.A.set(2.0);
    HeatSolver.VHC.set(3.0);
    HeatSolver.k.set(4.0);

    {
      ofstream output;
      output.open("Tbefore.txt");
      output << HeatSolver.T;
      output.close();
    }

    HeatSolver.stepForward(0.05);

    {
      ofstream output;
      output.open("Tafter.txt");
      output << HeatSolver.T;
      output.close();
    }

    {
      ofstream output;
      output.open("Tdiff.txt");
      output << Tinit;
      output.close();
    }
}
