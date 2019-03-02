#include "catch.hpp"

#include<iostream>
#include<fstream>

#include <LaserTherm/Emitter.hpp>
#include <Field.hpp>
#include <LaserTherm/HeatSources/_1D/BeersLaw.hpp>
#include <LaserTherm/Waveforms/ContinuousWave.hpp>
#include <LaserTherm/Waveforms/SinglePulse.hpp>
#include <LaserTherm/Waveforms/RegularPulseTrain.hpp>


TEST_CASE("Linear Absorption Emitter")
{

  SECTION("CW Exposure")
  {
  Emitter< HeatSources::_1D::BeersLaw<double>, Waveforms::ContinuousWave<double> > emitter(11);

  emitter.A.setCoordinateSystem( Uniform(0,2) );
  emitter.E0 = 2;
  emitter.mu_a.set(10);
  emitter.setStartTime(2);
  emitter.setExposureDuration(20);

  Field<double,1> A(11);
  A.set(0);

  emitter.addSourceTerm(A,0);

  CHECK(A(0)  == Approx(0));
  CHECK(A(5)  == Approx(0));
  CHECK(A(10) == Approx(0));

  emitter.addSourceTerm(A,2);

  CHECK(A(0)  == Approx(20));
  CHECK(A(5)  == Approx(20*exp(-10)));
  CHECK(A(10) == Approx(20*exp(-20)));

  emitter.addSourceTerm(A,20);

  CHECK(A(0)  == Approx(40));
  CHECK(A(5)  == Approx(40*exp(-10)));
  CHECK(A(10) == Approx(40*exp(-20)));

  emitter.addSourceTerm(A,30);

  CHECK(A(0)  == Approx(40));
  CHECK(A(5)  == Approx(40*exp(-10)));
  CHECK(A(10) == Approx(40*exp(-20)));
  }

  SECTION("Multi-Pulse Exposure")
  {

  Emitter< HeatSources::_1D::BeersLaw<double>, Waveforms::RegularPulseTrain<double> > emitter(11);

  emitter.A.setCoordinateSystem( Uniform(0,2) );
  emitter.E0 = 2;
  emitter.mu_a.set(10);
  emitter.setStartTime(2);
  emitter.setExposureDuration(20);
  emitter.setPulseDuration(1);
  emitter.setPulseRepititionFrequency(0.5);

  Field<double,1> A(11);
  A.set(0);

  emitter.addSourceTerm(A,0);

  CHECK(A(0)  == Approx(0));
  CHECK(A(5)  == Approx(0));
  CHECK(A(10) == Approx(0));

  emitter.addSourceTerm(A,2);

  CHECK(A(0)  == Approx(20));
  CHECK(A(5)  == Approx(20*exp(-10)));
  CHECK(A(10) == Approx(20*exp(-20)));

  emitter.addSourceTerm(A,3.1);

  CHECK(A(0)  == Approx(20));
  CHECK(A(5)  == Approx(20*exp(-10)));
  CHECK(A(10) == Approx(20*exp(-20)));

  emitter.addSourceTerm(A,20);

  CHECK(A(0)  == Approx(40));
  CHECK(A(5)  == Approx(40*exp(-10)));
  CHECK(A(10) == Approx(40*exp(-20)));

  emitter.addSourceTerm(A,21.1);

  CHECK(A(0)  == Approx(40));
  CHECK(A(5)  == Approx(40*exp(-10)));
  CHECK(A(10) == Approx(40*exp(-20)));

  emitter.addSourceTerm(A,30);

  CHECK(A(0)  == Approx(40));
  CHECK(A(5)  == Approx(40*exp(-10)));
  CHECK(A(10) == Approx(40*exp(-20)));

  }


}
