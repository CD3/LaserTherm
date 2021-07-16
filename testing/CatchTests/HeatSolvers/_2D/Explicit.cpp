#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

//#include <libField/HDF5.hpp>
#include <LaserTherm/HeatSolvers/_2D/Cylindrical/UniformExplicit.hpp>
#include <LaserTherm/HeatSolvers/_2D/Cylindrical/NonUniformExplicit.hpp>
using namespace std;

TEST_CASE("Explicit 2D Cylindrical Heat Solver","[heatsolver]")
{
    SECTION("Build UniformExplicit Heat Solver"){
      //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
      UniformExplicit<double> HeatSolver(20,10);
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

    SECTION("Build NonUniformExplicit Heat Solver"){
      //HeatSolvers::_2D::Cylindrical::Explicit<double> HeatSolver(10,20);
      NonUniformExplicit<double> HeatSolver(20,10);
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


TEST_CASE("UniformExplicit 2D Cylindrical Heat Solver Validation","[heatsolver][validation]")
{
  // see ./doc/writups/Validation/AnalyticalSolutions/AnalyticalSolutions.pdf
  // for a derivation of these tests
  int rN = 100;
  int zN = 200;
  double R = 2; // m
  double L = 4; // m
  double k = 0.5; // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18; // J/kg/K
  double dR = R / rN; // m
  double dL = L / zN; // m

  UniformExplicit<double> HeatSolver(zN-2,rN-2);
  Field<double, 2> Aplot(zN,rN);

  HeatSolver.T.setCoordinateSystem( Uniform(dL, L-dL), Uniform(dR, R-dR) );
  Aplot.setCoordinateSystem( Uniform(0.,L), Uniform(0., R) );


  HeatSolver.A.set(0.0);
  HeatSolver.VHC.set(rho*c);
  HeatSolver.k.set(k);
  HeatSolver.minZBC.f = 0;
  HeatSolver.maxZBC.f = 0;
  HeatSolver.maxRBC.f = 0;

  SECTION("Dirichlet BC")
  {
    double dt = 0.001;
    int Nt = 100;

    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 2.4048/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * sin(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
    };


    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });

    for(int i = 0; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }


    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");


    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }
  }

  SECTION("Nuemann BC"){
    // Neumann is --sad-- happy because his boundary condition test case is ---empty--- full --:(-- :)
    double dt = 0.001;
    int Nt = 100;

    HeatSolver.minZBC.type = BC::Type::HeatFlux;
    HeatSolver.maxZBC.type = BC::Type::HeatFlux;
    HeatSolver.maxRBC.type = BC::Type::HeatFlux;
    HeatSolver.minZBC.f = 0;
    HeatSolver.maxZBC.f = 0;
    HeatSolver.maxRBC.f = 0;

    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 3.8317/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * cos(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
    };

    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });
    for(int i = 0; i < Nt / 2; i++){
      HeatSolver.stepForward(dt);
    }
    for(int i = Nt / 2; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }

    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");


    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }
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

    for(int i = 0; i < Nt / 2; i++){
      HeatSolver.stepForward(dt);
    }



    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");

    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      //CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }

  }
}

TEST_CASE("NonUniformExplicit 2D Cylindrical Heat Solver Validation","[heatsolver][validation]")
{
  // see ./doc/writups/Validation/AnalyticalSolutions/AnalyticalSolutions.pdf
  // for a derivation of these tests
  int rN0 = 100;
  int zN0 = 200;
  int rN, zN;
  double r_s = 1.2;
  double z_s = 1.2;
  double R = 2; // m
  double L = 4; // m
  double k = 0.5; // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18; // J/kg/K
  double dR = R / rN0; // m
  double dL = L / zN0; // m
  // need to calculate n from finite geometric series formula
  rN = static_cast<int>(log( (1 - ( R * (1 - r_s) ) ) / dR ) / log(r_s));
  zN = static_cast<int>(log( (1 - ( L * (1 - z_s) ) ) / dL ) / log(z_s));

  NonUniformExplicit<double> HeatSolver(zN,rN);
  Field<double, 2> Aplot(zN,rN);

  HeatSolver.T.setCoordinateSystem( Geometric(0., dL, z_s), Geometric(0., dR, r_s) );
  Aplot.setCoordinateSystem( Geometric(0., dL, z_s), Geometric(0., dR, r_s) );


  HeatSolver.A.set(0.0);
  HeatSolver.VHC.set(rho*c);
  HeatSolver.k.set(k);
  HeatSolver.minZBC.f = 0;
  HeatSolver.maxZBC.f = 0;
  HeatSolver.maxRBC.f = 0;

  SECTION("NonUniform Dirichlet BC")
  {
    double dt = 0.001;
    int Nt = 100;

    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 2.4048/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * sin(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
    };

    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });

    for(int i = 0; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }

    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");

    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first, temp.second)).epsilon(0.01));
    }
  }

  SECTION("Nuemann BC"){
    // Neumann is --sad-- happy because his boundary condition test case is ---empty--- full --:(-- :)
    double dt = 0.001;
    int Nt = 100;

    HeatSolver.minZBC.type = BC::Type::HeatFlux;
    HeatSolver.maxZBC.type = BC::Type::HeatFlux;
    HeatSolver.maxRBC.type = BC::Type::HeatFlux;
    HeatSolver.minZBC.f = 0;
    HeatSolver.maxZBC.f = 0;
    HeatSolver.maxRBC.f = 0;

    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 3.8317/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * cos(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
    };

    HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    Aplot.set_f([&](auto x) { return solution(x[0],x[1],Nt*dt); });
    for(int i = 0; i < Nt / 2; i++){
      HeatSolver.stepForward(dt);
    }
    for(int i = Nt / 2; i < Nt; i++){
      HeatSolver.stepForward(dt);
    }
    {
      ofstream output;
      output.open("T1.txt");
      output << HeatSolver.T;
      output.close();
    }
    {
      ofstream output;
      output.open("T2.txt");
      output << Aplot;
      output.close();
    }

    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");


    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first, temp.second)).epsilon(0.01));
    }
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

    for(int i = 0; i < Nt / 2; i++){
      HeatSolver.stepForward(dt);
    }



    vector<string> Name;
    Name.push_back("center");
    Name.push_back("rmax");
    Name.push_back("zmax");
    Name.push_back("zmin");
    Name.push_back("r=0");

    vector<std::pair<int, int>> Points;
    Points.push_back(std::make_pair<int, int>(zN/2, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, rN-5));
    Points.push_back(std::make_pair<int, int>(zN-5, rN/2));
    Points.push_back(std::make_pair<int, int>(3, rN/2));
    Points.push_back(std::make_pair<int, int>(zN/2, 3));

    for(int i = 0; i < Points.size(); i++){
      std::pair<int, int> temp = Points[i];
      //std::cout << "Testing " << Name[i] << "\n";
      //CHECK( HeatSolver.T(temp.first, temp.second) == Approx(Aplot(temp.first+1, temp.second+1)).epsilon(0.01));
    }

  }
}


TEST_CASE("Auto Adaptive Time Step","[heatsolver][validation]"){
  int rN0 = 100;
  int zN0 = 200;
  int rN, zN;
  double r_s = 1.18;
  double z_s = 1.18;
  double R = 2; // m
  double L = 4; // m
  double k = 0.5; // W/m/K
  double rho = 1000; // kg/m^3
  double c = 4.18; // J/kg/K
  double dR = R / rN0; // m
  double dL = L / zN0; // m
  // need to calculate n from finite geometric series formula
  rN = static_cast<int>(log( (1 - ( R * (1 - r_s) ) / dR ) ) / log(r_s)) + 2;
  zN = static_cast<int>(log( (1 - ( L * (1 - z_s) ) / dL ) ) / log(z_s)) + 2;

  NonUniformExplicit<double> NU_HeatSolver(zN,rN);
  UniformExplicit<double> U_HeatSolver(200,100);
  Field<double, 2> NU_Analytic(zN,rN);
  Field<double, 2> U_Analytic(zN,rN);

  NU_HeatSolver.T.setCoordinateSystem( Geometric(0., dL, z_s), Geometric(0., dR, r_s) );
  U_HeatSolver.T.setCoordinateSystem( Uniform(0., L), Uniform(0., R) );
  NU_Analytic.setCoordinateSystem( Geometric(0., dL, z_s), Geometric(0., dR, r_s) );
  U_Analytic.setCoordinateSystem( Uniform(0., L), Uniform(0., R) );


  NU_HeatSolver.A.set(0.0);
  NU_HeatSolver.VHC.set(rho*c);
  NU_HeatSolver.k.set(k);
  NU_HeatSolver.minZBC.f = 0;
  NU_HeatSolver.maxZBC.f = 0;
  NU_HeatSolver.maxRBC.f = 0;

  U_HeatSolver.A.set(0.0);
  U_HeatSolver.VHC.set(rho*c);
  U_HeatSolver.k.set(k);
  U_HeatSolver.minZBC.f = 0;
  U_HeatSolver.maxZBC.f = 0;
  U_HeatSolver.maxRBC.f = 0;

  SECTION("Algorithm Convergence")
  {
    double NUdt = 100.;
    double Udt = 100.;
    double t_f = 10.;

    auto solution = [&](double z, double r, double t)
    {
      double lambda_r = 2.4048/R;
      double lambda_z = M_PI/L;
      double alpha = k/rho/c*( lambda_z*lambda_z + lambda_r*lambda_r);

      return exp(-alpha*t) * sin(lambda_z*z) * std::cyl_bessel_j(0,lambda_r*r);
    };

    NU_HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    U_HeatSolver.T.set_f([&](auto x) { return solution(x[0],x[1],0); });
    NU_Analytic.set_f([&](auto x) { return solution(x[0],x[1],t_f); });
    U_Analytic.set_f([&](auto x) { return solution(x[0],x[1],t_f); });

    NUdt = NU_HeatSolver.findTimeStep(NUdt, 0.0, 0.01);
    std::cout << "NonUniform Time Step: " << NUdt << "\n";

    Udt = U_HeatSolver.findTimeStep(Udt, 0.0, 0.01);
    std::cout << "Uniform Time Step: " << Udt << "\n";

    double U_err = 0;
    double NU_err = 0;
    for(int i = 0; i < U_Analytic.size(0); i++){
      for(int j = 0; j < U_Analytic.size(1); j++){
        if(U_HeatSolver.T[i][j] + U_Analytic[i][j] == 0){
          if(U_HeatSolver.T[i][j] + U_Analytic[i][j] == U_HeatSolver.T[i][j] - U_Analytic[i][j]){
            continue;
          } else{
            U_err += U_HeatSolver.T[i][j];
          }
        } else{
        U_err += (U_HeatSolver.T[i][j] - U_Analytic[i][j]) / (U_HeatSolver.T[i][j] + U_Analytic[i][j]); 
        }
      }
    }
    U_err /= U_Analytic.size(0) * U_Analytic.size(1);

    for(int i = 0; i < NU_Analytic.size(0); i++){
      for(int j = 0; j < NU_Analytic.size(1); j++){
        if(NU_HeatSolver.T[i][j] + NU_Analytic[i][j] == 0){
          if(NU_HeatSolver.T[i][j] + NU_Analytic[i][j] == NU_HeatSolver.T[i][j] - NU_Analytic[i][j]){
            continue;
          } else{
            NU_err += NU_HeatSolver.T[i][j];
          }
        } else{
        NU_err += (NU_HeatSolver.T[i][j] - NU_Analytic[i][j]) / (NU_HeatSolver.T[i][j] + NU_Analytic[i][j]); 
        }
      }
    }
    NU_err /= NU_Analytic.size(0) * NU_Analytic.size(1);

    std::cout << "Uniform Error:" << U_err << "\n";
    std::cout << "NonUniform Error:" << NU_err << "\n";

    //moveForward in time using an adaptive timestep
    NU_HeatSolver.moveForward(10., 0.0, 0.01);
    {
      ofstream output;
      output.open("NU_HeatSolACTS.txt");
      output << NU_HeatSolver.T;
      output.close();
    }
    U_HeatSolver.moveForward(10., 0.0, 0.01);
    {
      ofstream output;
      output.open("U_HeatSolACTS.txt");
      output << U_HeatSolver.T;
      output.close();
    }

  }
}
