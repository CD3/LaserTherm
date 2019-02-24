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
          BC.type = HeatSolvers::BoundaryConditions::Type::Temperature;
          BC.f = 10;
        }
    );

    HeatSolver.sig_askMaxBoundaryCondition.connect(
        [](auto& BC, const auto& T, const auto& t)
        {
          BC.type = HeatSolvers::BoundaryConditions::Type::HeatFlux;
          BC.f = 20;
          BC.dfdT = 2;
        }
    );

    CHECK( HeatSolver.minBC.type == HeatSolvers::BoundaryConditions::Type::None );
    CHECK( HeatSolver.maxBC.type == HeatSolvers::BoundaryConditions::Type::None );

    HeatSolver.sig_askMinBoundaryCondition( HeatSolver.minBC, 0, 0 );

    CHECK( HeatSolver.minBC.type == HeatSolvers::BoundaryConditions::Type::Temperature);
    CHECK( HeatSolver.maxBC.type == HeatSolvers::BoundaryConditions::Type::None );
    CHECK( HeatSolver.minBC.f == 10 );

    HeatSolver.sig_askMaxBoundaryCondition( HeatSolver.maxBC, 0, 0 );

    CHECK( HeatSolver.minBC.type == HeatSolvers::BoundaryConditions::Type::Temperature);
    CHECK( HeatSolver.maxBC.type == HeatSolvers::BoundaryConditions::Type::HeatFlux);
    CHECK( HeatSolver.minBC.f == 10 );
    CHECK( HeatSolver.maxBC.f == 20 );
    CHECK( HeatSolver.maxBC.dfdT == 2 );






  }
}

TEST_CASE("1D Cartesian Heat Solver Validation", "[heatsolver][validation]")
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
        HeatSolver.minBC.type = HeatSolvers::BoundaryConditions::Type::HeatFlux;
        HeatSolver.maxBC.type = HeatSolvers::BoundaryConditions::Type::HeatFlux;

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
        HeatSolver.minBC.type = HeatSolvers::BoundaryConditions::Type::HeatFlux;
        HeatSolver.maxBC.type = HeatSolvers::BoundaryConditions::Type::HeatFlux;

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

  SECTION("Boundary Conditions")
  {

      SECTION("Constant Temperature Boundaries")
      {
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> min(-1);
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> max(1);

        min.setBoundaryCondition( HeatSolver.minBC );
        max.setBoundaryCondition( HeatSolver.maxBC );

        auto sol = [&](int i){ return (2./L)*(xmin + i*dx) - 1; };

        // set the initial temperature distribution to something close to the solution so we
        // don't have to wait so long.
        HeatSolver.T.set_f(
            [&](auto i, auto cs)
            {
             return sol(i[0])+0.1;
            }
            );

        for(int i = 0; i < 10000; ++i)
          HeatSolver.stepForward( HeatSolver.calcMaxTimeStep() );

        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.01));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.01));
      }

      SECTION("Constant Heat Flux at Max Boundary")
      {
        double Q = 0.01;
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> min(0);
        HeatSolvers::BoundaryConditions::ConstantHeatFlux<double> max(-Q);

        min.setBoundaryCondition( HeatSolver.minBC );
        max.setBoundaryCondition( HeatSolver.maxBC );

        auto sol = [&](int i){ return -(Q/k)*(xmin + i*dx); };

        HeatSolver.T.set_f(
            [&](auto i, auto cs)
            {
             return sol(i[0])+0.1;
            }
            );

        for(int i = 0; i < 10000; ++i)
          HeatSolver.stepForward( HeatSolver.calcMaxTimeStep() );

        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.05));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.01));
      }

      SECTION("Constant Heat Flux at Min Boundary")
      {
        double Q = 0.01;
        HeatSolvers::BoundaryConditions::ConstantHeatFlux<double> min(-Q);
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> max(0);

        min.setBoundaryCondition( HeatSolver.minBC );
        max.setBoundaryCondition( HeatSolver.maxBC );


        auto sol = [&](int i){ return (Q/k)*(xmin + i*dx - L); };

        HeatSolver.T.set_f(
            [&](auto i, auto cs)
            {
             return sol(i[0])+0.1;
            }
            );

        for(int i = 0; i < 10000; ++i)
          HeatSolver.stepForward( HeatSolver.calcMaxTimeStep() );

        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.01));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.05));
      }

      SECTION("Convective Heat Flux at Max Boundary")
      {
        double Tinf = 10;
        double he = 0.01;
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> min(0);
        HeatSolvers::BoundaryConditions::Convective<double> max(Tinf,he);


        min.setBoundaryCondition( HeatSolver.minBC );
        HeatSolver.sig_askMaxBoundaryCondition.connect( [&max](auto& BC, const auto& T, const auto& t){
            max.setBoundaryCondition(BC,T,t);
            } );


        auto sol = [&](int i){ return (he*Tinf/(k + he*L))*(xmin + i*dx); };

        HeatSolver.T.set_f(
            [&](auto i, auto cs)
            {
             return sol(i[0])+0.1;
            }
            );

        for(int i = 0; i < 10000; ++i)
        {
          HeatSolver.stepForward( HeatSolver.calcMaxTimeStep() );
        }

        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.05));
        CHECK( HeatSolver.T(7*N/8) == Approx( sol(7*N/8) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N-1) == Approx( sol(  N-1) ).epsilon(0.01));
      }

      SECTION("Convective Heat Flux at Min Boundary")
      {
        double Tinf = 10;
        double he = 0.01;
        HeatSolvers::BoundaryConditions::ConstantTemperature<double> max(0);
        HeatSolvers::BoundaryConditions::Convective<double> min(Tinf,he);


        max.setBoundaryCondition( HeatSolver.maxBC );
        HeatSolver.sig_askMinBoundaryCondition.connect( [&min](auto& BC, const auto& T, const auto& t){
            min.setBoundaryCondition(BC,T,t);
            } );


        auto sol = [&](int i){ return (-he*Tinf/(k + he*L))*(xmin + i*dx - L); };

        HeatSolver.T.set_f(
            [&](auto i, auto cs)
            {
             return sol(i[0])+0.1;
            }
            );

        std::ofstream out("T.txt");
        out << HeatSolver.T;
        for(int i = 0; i < 10000; ++i)
        {
          HeatSolver.stepForward( HeatSolver.calcMaxTimeStep() );
          out << "\n" << HeatSolver.T;
        }
        CHECK( HeatSolver.T(    0) == Approx( sol(    0) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/2) == Approx( sol(  N/2) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/4) == Approx( sol(  N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(3*N/4) == Approx( sol(3*N/4) ).epsilon(0.01));
        CHECK( HeatSolver.T(  N/8) == Approx( sol(  N/8) ).epsilon(0.05));
      }


    }

  }

}
