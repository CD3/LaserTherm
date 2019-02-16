#include "catch.hpp"

#include <iostream>
#include <fstream>
#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>

template<typename T>
void nop( T& a )
{
}

TEST_CASE("Construction and Setup", "[heatsolver]")
{

  HeatSolvers::_1D::CrankNicholson<double> HeatSolver(11);
  HeatSolver.T.setCoordinateSystem( Uniform(-2,2) );

  HeatSolver.T.set(1.0);
  HeatSolver.A.set(2.0);
  HeatSolver.VHC.set(3.0);
  HeatSolver.k.set(4.0);

  CHECK( HeatSolver.T.getCoord(0) == Approx(-2) );
  CHECK( HeatSolver.T.getCoord(1) == Approx(-1.6) );
  CHECK( HeatSolver.T.getCoord(9) == Approx(1.6) );
  CHECK( HeatSolver.T.getCoord(10) == Approx(2) );


  CHECK( HeatSolver.A.getCoord(0) == Approx(-2) );
  CHECK( HeatSolver.A.getCoord(1) == Approx(-1.6) );
  CHECK( HeatSolver.A.getCoord(9) == Approx(1.6) );
  CHECK( HeatSolver.A.getCoord(10) == Approx(2) );


  CHECK( HeatSolver.VHC.getCoord(0) == Approx(-2) );
  CHECK( HeatSolver.VHC.getCoord(1) == Approx(-1.6) );
  CHECK( HeatSolver.VHC.getCoord(9) == Approx(1.6) );
  CHECK( HeatSolver.VHC.getCoord(10) == Approx(2) );

  CHECK( HeatSolver.k.getCoord(0) == Approx(-2) );
  CHECK( HeatSolver.k.getCoord(1) == Approx(-1.6) );
  CHECK( HeatSolver.k.getCoord(9) == Approx(1.6) );
  CHECK( HeatSolver.k.getCoord(10) == Approx(2) );

  CHECK( HeatSolver.T(0) == 1.0 );
  CHECK( HeatSolver.A(0) == 2.0 );
  CHECK( HeatSolver.VHC(0) == 3.0 );
  CHECK( HeatSolver.k(0) == 4.0 );

  CHECK( HeatSolver.T(10) == 1.0 );
  CHECK( HeatSolver.A(10) == 2.0 );
  CHECK( HeatSolver.VHC(10) == 3.0 );
  CHECK( HeatSolver.k(10) == 4.0 );


}

TEST_CASE("1D Cartesian Heat Solver Validation", "[heatsolver,validation]")
{

  // for thermally homogeneous material of length L and sink boundary conditions, the eigenfunctions of the heat conduction operator are sin functions.
  // these eigen functions then decay exponentially and their decay rate is given by their frequency (and the thermal properties).
  // so, if we load an initial temperature profile that is a sin wave, the shape of the profile will not change, and the whole thing will
  // exponentially decay.
  //
  //
  //
  //   |---------------------------------------------|
  //   0                                             L
  //
  //
  //

  int N = 600;

  double L  = 3.0;
  double xmin = 0;
  double xmax = L;
  double dx = (xmax-xmin)/(N-1);

  double rho = 1.0;
  double c = 4.1868;
  double k = 0.00628;

  int    Nt = 1000;
  double dt = 0.1;
  double Tmax = Nt*dt;


  int m = 2;
  double alpha = pow( m * M_PI / L, 2) * k / (rho*c); // decay rate for the mode


  SECTION("Crank-Nicholson Solver")
  {

    HeatSolvers::_1D::CrankNicholson<double> HeatSolver(N);

    SECTION("Sink Boundary Conditions")
    {

      SECTION("Direct Field Access")
      {
        HeatSolver.T.setCoordinateSystem( Uniform(xmin,xmax) );
        HeatSolver.VHC.set( rho*c );

        HeatSolver.k.set( k );

        CHECK( HeatSolver.calcMaxTimeStep() == Approx(dx * rho * c / 2 / k) );

        HeatSolver.T.set_f( [&](auto i, auto cs)
        {
          auto x = cs->getCoord(i[0]);
          return sin( m * M_PI * x / L );
        });

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }

      }

      SECTION("Callback Field Access")
      {
        HeatSolver.T.setCoordinateSystem( Uniform(xmin,xmax) );

        HeatSolver.sig_askInitialTemperature.connect( [&](Field<double,1>& T)
            {
              size_t N = T.size(0);
              for(int i = 0; i < N; ++i)
              {
                double x = HeatSolver.T.getCoord(i);
                T(i) = sin( m * M_PI * x / L );
              }
            });

        HeatSolver.sig_askVolumetricHeatCapacity.connect( [&](Field<double,1>& VHC)
            {
              VHC.set( rho*c );
            });

        HeatSolver.sig_askConductivity.connect( [&](Field<double,1>& kk)
            {
              kk.set( k );
            });


        HeatSolver.initialize();
        for(int i = 0; i < N; ++i)
        {
          CHECK( HeatSolver.VHC(i) == Approx(rho*c) );
          CHECK( HeatSolver.k(i) == Approx(k) );
          CHECK( HeatSolver.T(i) == Approx( sin(m*M_PI*(xmin+i*dx)/L) ) );
        }

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }

      }

      CHECK( HeatSolver.T(  N/2) == Approx( exp( - alpha * Tmax) * sin( m * M_PI * (xmin + dx*  N/2) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(  N/8) == Approx( exp( - alpha * Tmax) * sin( m * M_PI * (xmin + dx*  N/8) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(7*N/8) == Approx( exp( - alpha * Tmax) * sin( m * M_PI * (xmin + dx*7*N/8) / L )  ).epsilon(0.01));
    }

    SECTION("Insulating Boundary Conditions")
    {
      HeatSolver.T.setCoordinateSystem( Uniform(xmin,xmax) );
      HeatSolver.VHC.set( rho*c );
      HeatSolver.k.set( k );
      HeatSolver.minBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;
      HeatSolver.maxBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;

      for(int i = 0; i < N; ++i)
      {
        double x = HeatSolver.T.getCoord(i);
        HeatSolver.T(i) = cos( m * M_PI * x / L );
      }

      std::ofstream out("T.txt");
      out << HeatSolver.T;
      for(int i = 0; i < Nt; ++i)
      {
        HeatSolver.stepForward(dt);
        out << "\n" << HeatSolver.T;
      }
      out.close();


      CHECK( HeatSolver.T(    0) == Approx( exp( - alpha * Tmax) * cos( m * M_PI * (xmin + dx*    0) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(  N/2) == Approx( exp( - alpha * Tmax) * cos( m * M_PI * (xmin + dx*  N/2) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(  N/8) == Approx( exp( - alpha * Tmax) * cos( m * M_PI * (xmin + dx*  N/8) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(7*N/8) == Approx( exp( - alpha * Tmax) * cos( m * M_PI * (xmin + dx*7*N/8) / L )  ).epsilon(0.01));
      CHECK( HeatSolver.T(  N-1) == Approx( exp( - alpha * Tmax) * cos( m * M_PI * (xmax           ) / L )  ).epsilon(0.01));
    }

  }




}
