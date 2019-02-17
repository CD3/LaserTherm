#include "catch.hpp"

#include <iostream>
#include <fstream>
#include <LaserTherm/HeatSolvers/_1D/CrankNicholson.hpp>

template<typename T>
void nop( T& a )
{
}

TEST_CASE("Heat Solver Construction and Setup", "[heatsolver]")
{

  SECTION("Crank-Nicholson")
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


}

TEST_CASE("Heat Solver Signals", "[heatsolver]")
{
  SECTION("Crank-Nicholson")
  {
    HeatSolvers::_1D::CrankNicholson<double> HeatSolver(10);

    HeatSolver.sig_askSourceTerm.connect(
        [](auto& f, const auto& t)
        {
          f.set(2);

          if( t > 10 )
            f.set(10);

          if( t > 20 )
            f.set(20);

          if( t < 0.5 )
            f.set(1);
        }
    );

    HeatSolver.sig_askSourceTerm( HeatSolver.A, 11 );
    CHECK( HeatSolver.A(0) == 10 );
    CHECK( HeatSolver.A(1) == 10 );
    CHECK( HeatSolver.A(8) == 10 );
    CHECK( HeatSolver.A(9) == 10 );

    HeatSolver.sig_askSourceTerm( HeatSolver.A, 21 );
    CHECK( HeatSolver.A(0) == 20 );
    CHECK( HeatSolver.A(1) == 20 );
    CHECK( HeatSolver.A(8) == 20 );
    CHECK( HeatSolver.A(9) == 20 );

    HeatSolver.stepForward(0.1);
    CHECK( HeatSolver.A(0) == 1 );
    CHECK( HeatSolver.A(1) == 1 );
    CHECK( HeatSolver.A(8) == 1 );
    CHECK( HeatSolver.A(9) == 1 );

    HeatSolver.stepForward(1);
    CHECK( HeatSolver.A(0) == 2 );
    CHECK( HeatSolver.A(1) == 2 );
    CHECK( HeatSolver.A(8) == 2 );
    CHECK( HeatSolver.A(9) == 2 );


    HeatSolver.sig_askMinBoundaryCondition.connect(
        [](auto& BC, const auto& T, const auto& t)
        {
          BC.type = HeatSolvers::_1D::BoundaryConditions::Type::Dirichlet;
          BC.f = 10;
        }
    );

    HeatSolver.sig_askMaxBoundaryCondition.connect(
        [](auto& BC, const auto& T, const auto& t)
        {
          BC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;
          BC.f = 20;
          BC.dfdT = 2;
        }
    );

    CHECK( HeatSolver.minBC.type == HeatSolvers::_1D::BoundaryConditions::Type::None );
    CHECK( HeatSolver.maxBC.type == HeatSolvers::_1D::BoundaryConditions::Type::None );

    HeatSolver.sig_askMinBoundaryCondition( HeatSolver.minBC, 0, 0 );

    CHECK( HeatSolver.minBC.type == HeatSolvers::_1D::BoundaryConditions::Type::Dirichlet);
    CHECK( HeatSolver.maxBC.type == HeatSolvers::_1D::BoundaryConditions::Type::None );
    CHECK( HeatSolver.minBC.f == 10 );

    HeatSolver.sig_askMaxBoundaryCondition( HeatSolver.maxBC, 0, 0 );

    CHECK( HeatSolver.minBC.type == HeatSolvers::_1D::BoundaryConditions::Type::Dirichlet);
    CHECK( HeatSolver.maxBC.type == HeatSolvers::_1D::BoundaryConditions::Type::Neumann);
    CHECK( HeatSolver.minBC.f == 10 );
    CHECK( HeatSolver.maxBC.f == 20 );
    CHECK( HeatSolver.maxBC.dfdT == 2 );






  }
}

TEST_CASE("1D Cartesian Heat Solver Validation", "[heatsolver,validation]")
{
  // see ./doc/writups/Validation/AnalyticalSolutions/AnalyticalSolutions.pdf for a derivation
  // of these tests

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
  double lambda = m * M_PI / L;
  double alpha = k * pow( lambda, 2) / (rho*c); // decay rate for the mode
  double beta = 5;


  SECTION("Crank-Nicholson Solver")
  {

    HeatSolvers::_1D::CrankNicholson<double> HeatSolver(N);

    HeatSolver.T.setCoordinateSystem( Uniform(xmin,xmax) );
    HeatSolver.VHC.set( rho*c );
    HeatSolver.k.set( k );

    CHECK( HeatSolver.calcMaxTimeStep() == Approx(dx * rho * c / 2 / k) );

    SECTION("Without Source")
    {

      SECTION("Sink Boundary Conditions")
      {

        HeatSolver.T.set_f( [&](auto i, auto cs)
        {
          auto x = cs->getCoord(i[0]);
          return sin( lambda * x );
        });

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }

        auto sol = [&](int i){ return exp( - alpha * Tmax) * sin( lambda * (xmin + dx*i) );};
        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.01));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.01));
      }

      SECTION("Insulating Boundary Conditions")
      {
        HeatSolver.minBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;
        HeatSolver.maxBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;

        for(int i = 0; i < N; ++i)
        {
          double x = HeatSolver.T.getCoord(i);
          HeatSolver.T(i) = cos( lambda * x );
        }

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }

        auto sol = [&](int i){ return exp( - alpha * Tmax) * cos( lambda * (xmin + dx*i) );};
        CHECK( HeatSolver.T(    0) == Approx( sol(    0) ).epsilon(0.001));
        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.001));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.001));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.001));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.001));
        CHECK( HeatSolver.T(  N-1) == Approx( sol(  N-1) ).epsilon(0.001));

      }

    }



  SECTION("With Source")
  {
      SECTION("Sink Boundary Conditions")
      {

        HeatSolver.T.set(0.0);

        HeatSolver.A.set_f( [&](auto i, auto cs)
        {
          auto x = cs->getCoord(i[0]);
          return beta*sin( lambda * x );
        });

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }

        auto sol = [&](int i, int n){ return (beta/k/lambda/lambda)*(1 - exp( -alpha*dt*n)) * sin( lambda * (xmin + dx*i) );};
        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8,Nt) ).epsilon(0.01));
      }

      SECTION("Insulating Boundary Conditions")
      {
        HeatSolver.minBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;
        HeatSolver.maxBC.type = HeatSolvers::_1D::BoundaryConditions::Type::Neumann;

        HeatSolver.A.set_f( [&](auto i, auto cs)
        {
          auto x = cs->getCoord(i[0]);
          return beta*cos( lambda * x );
        });

        for(int i = 0; i < Nt; ++i)
        {
          HeatSolver.stepForward(dt);
        }


        auto sol = [&](int i, int n){ return (beta/k/lambda/lambda)*(1 - exp( -alpha*dt*n)) * cos( lambda * (xmin + dx*i) );};
        CHECK( HeatSolver.T(    0) == Approx( sol(    0,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8,Nt) ).epsilon(0.01));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8,Nt) ).epsilon(0.01));

      }




    }

  }

}
