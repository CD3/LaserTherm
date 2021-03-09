#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

#include <libField/HDF5.hpp>
#include <LaserTherm/HeatSolvers/_2D/Cylindrical/Explicit.hpp>
using namespace std;

TEST_CASE("Explicit 2D Cylindrical Heat Solver","[heatsolver]")
{
    SECTION("Sandbox testing area"){
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
    SECTION("Build Explicit Heat Solver"){
      // Create and set variables for Explicit. Checking constructors and accessing members
    }
    SECTION("Step Forward"){
      // Check step Forward runs successfully, keeps a constant field constant,
      // and the center point of a gaussian GOES DOWN!
    }
    SECTION("Boundary Conditions"){
      // Check boundary conditions hold up somehow... what's a good test?
    }

}

TEST_CASE("Analytic Comparisons"){
  SECTION("Rod with end at constant temp"){
  
  }
  SECTION("Bessel Function"){
  
  }
}

TEST_CASE("Explicit 2D Cylindrical Heat Solver Validation","[heatsolver][validation]")
{
  // see ./doc/writups/Validation/AnalyticalSolutions/AnalyticalSolutions.pdf
  // for a derivation of these tests
  Explicit<double> HeatSolver(100,200);
  double R = 2;
  double L = 4;
  HeatSolver.T.setCoordinateSystem( Uniform(0.,R), Uniform(0., L) );

  double k = 0.5; // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18; // J/kg/K
  HeatSolver.A.set(0.0);
  HeatSolver.VHC.set(rho*c);
  HeatSolver.k.set(k);

  SECTION("Dirichlet BC")
  {
    auto solution = [&](double r, double z, double t)
    {
      double lambda_r = 2.4048/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * sin(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
      
    };

    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });

    /* hdf5write("T_initial.h5", HeatSolver.T); */

    double dt = 0.01;
    int Nt = 100;
    for(int i = 0; i < Nt; i++){
      HeatSolver.stepForward(0.01);
    }

    CHECK( HeatSolver.T(50,1) == Approx( solution(R/2,L/2,Nt*dt) ) );


  }
}
