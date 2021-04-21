#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

//#include <libField/HDF5.hpp>
#include <LaserTherm/HeatSolvers/_2D/Cylindrical/Explicit.hpp>
using namespace std;

TEST_CASE("Explicit 2D Cylindrical Heat Solver","[heatsolver]")
{
    SECTION("Build Explicit Heat Solver"){
      //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
      Explicit<double> HeatSolver(20,10);
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

      HeatSolver.stepForward(0.001);
    }
}


TEST_CASE("Explicit 2D Cylindrical Heat Solver Validation","[heatsolver][validation]")
{
  // see ./doc/writups/Validation/AnalyticalSolutions/AnalyticalSolutions.pdf
  // for a derivation of these tests
  double R = 2;
  double L = 4;
  int rN = 141;
  int zN = 282;
  double dR = R / rN;
  double dL = L / zN;
  Explicit<double> HeatSolver(zN-2,rN-2);
  Field<double, 2> Aplot(zN,rN);
  HeatSolver.T.setCoordinateSystem( Uniform(dL, L-dL), Uniform(dR, R-dR) );
  Aplot.setCoordinateSystem( Uniform(0.,L), Uniform(0., R) );


  double k = 0.5; // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18; // J/kg/K
  HeatSolver.A.set(0.0);
  HeatSolver.VHC.set(rho*c);
  HeatSolver.k.set(k);
  HeatSolver.minZBC.f = 0;
  HeatSolver.maxZBC.f = 0;
  HeatSolver.maxRBC.f = 0;

  SECTION("Dirichlet BC")
  {
    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 2.4048/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * sin(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);

    };

    double dt = 0.001;
    int Nt = 100;
    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });
    for(int i = 0; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }

    {
      ofstream output;
      output.open("NumericalSol.txt");
      output << HeatSolver.T;
      output.close();
    }
    {
      ofstream output;
      output.open("AnalyticalSol.txt");
      output << Aplot;
      output.close();
    }

    /* hdf5write("T_initial.h5", HeatSolver.T); */

    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");


    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(100, 50));
    Points.push_back(std::make_pair<int, int>(100, 98-1));
    Points.push_back(std::make_pair<int, int>(198-1, 50));
    Points.push_back(std::make_pair<int, int>(1, 50));
    Points.push_back(std::make_pair<int, int>(100, 1));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      std::cout << "Testing " << Name[i] << "\n";
      CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }
  }

  SECTION("Nuemann BC"){
    // Neumann is sad because his boundary condition test case is empty :(
  }

  SECTION("Green's Function"){
    // Remember: reset values of HeatSolver(/aplot)
    auto solution = [&](double z, double r, double t)
    {
      return 0;

    };

    double dt = 0.001;
    int Nt = 100;

    HeatSolver.T.set(0);
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });

    HeatSolver.A[zN / 2][rN / 2] = 1;

    for(int i = 0; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }

    {
      ofstream output;
      output.open("NumericalSol.txt");
      output << HeatSolver.T;
      output.close();
    }
    {
      ofstream output;
      output.open("AnalyticalSol.txt");
      output << Aplot;
      output.close();
    }

    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");

    CHECK( HeatSolver.T(zN / 2, rN / 2) == Approx(1));

    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(100, 50));
    Points.push_back(std::make_pair<int, int>(100, 98-1));
    Points.push_back(std::make_pair<int, int>(198-1, 50));
    Points.push_back(std::make_pair<int, int>(1, 50));
    Points.push_back(std::make_pair<int, int>(100, 1));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      std::cout << "Testing " << Name[i] << "\n";
      //CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }

  }
}
