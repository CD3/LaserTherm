#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>

#include <LaserTherm/Waveforms/Mixins.hpp>

TEST_CASE("Waveform Mixins") {
  SECTION("CW Waveform Storage") {
    Waveforms::AddExposureDuration<
        Waveforms::AddStartTime<Waveforms::Base<double>>>
        exp;

    exp.setExposureDuration(10);
    exp.setStartTime(2);

    CHECK(exp.getExposureDuration() == Approx(10));
    CHECK(exp.getStartTime() == Approx(2));
  }

  SECTION("CW Waveform isOn Check") {
    SECTION("with start time and exposure duration") {
      Waveforms::AddIsOnCalc<Waveforms::AddExposureDuration<
          Waveforms::AddStartTime<Waveforms::Base<double>>>>
          exp;

      exp.setExposureDuration(10);
      exp.setStartTime(2);

      CHECK(!exp.isOn(0.0));
      CHECK(!exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(exp.isOn(11.0));
      CHECK(exp.isOn(12.0));
      CHECK(!exp.isOn(13.0));
      CHECK(!exp.isOn(20.));
    }

    SECTION("with exposure duration but no start time") {
      Waveforms::AddIsOnCalc<
          Waveforms::AddExposureDuration<Waveforms::Base<double>>>
          exp;

      exp.setExposureDuration(10);

      CHECK(exp.isOn(0.0));
      CHECK(exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(!exp.isOn(11.0));
      CHECK(!exp.isOn(12.0));
      CHECK(!exp.isOn(13.0));
      CHECK(!exp.isOn(20.));
    }

    SECTION("with start time but no exposure duration") {
      Waveforms::AddIsOnCalc<Waveforms::AddStartTime<Waveforms::Base<double>>>
          exp;

      exp.setStartTime(2);

      CHECK(!exp.isOn(0.0));
      CHECK(!exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(exp.isOn(11.0));
      CHECK(exp.isOn(12.0));
      CHECK(exp.isOn(13.0));
      CHECK(exp.isOn(20.));
    }

    SECTION("with no start time or exposure duration") {
      Waveforms::AddIsOnCalc<Waveforms::Base<double>> exp;

      CHECK(exp.isOn(0.0));
      CHECK(exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(exp.isOn(11.0));
      CHECK(exp.isOn(12.0));
      CHECK(exp.isOn(13.0));
      CHECK(exp.isOn(20.));
    }
  }

  SECTION("Regular Pulse Train Storage") {
    Waveforms::AddExposureDuration<
        Waveforms::AddStartTime<Waveforms::AddPulseDuration<
            Waveforms::AddPulseRepititionFrequency<Waveforms::Base<double>>>>>
        exp;

    exp.setExposureDuration(10);
    exp.setStartTime(2);
    exp.setPulseDuration(0.01);
    exp.setPulseRepititionFrequency(10);

    CHECK(exp.getExposureDuration() == Approx(10));
    CHECK(exp.getStartTime() == Approx(2));
    CHECK(exp.getPulseDuration() == Approx(0.01));
    CHECK(exp.getPulseRepititionFrequency() == Approx(10));
  }

  SECTION("Regular Pulse Train Waveform isOn Check") {
    SECTION("with start time, exposure duration, pulse duration, and PRF") {
      Waveforms::AddIsOnCalc<Waveforms::AddExposureDuration<
          Waveforms::AddStartTime<Waveforms::AddPulseDuration<
              Waveforms::AddPulseRepititionFrequency<
                  Waveforms::Base<double>>>>>>
          exp;

      exp.setExposureDuration(10);
      exp.setStartTime(2);
      exp.setPulseDuration(0.01);
      exp.setPulseRepititionFrequency(10);

      CHECK(!exp.isOn(0.0));
      CHECK(!exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(exp.isOn(11.0));
      CHECK(exp.isOn(12.0));
      CHECK(!exp.isOn(13.0));
      CHECK(!exp.isOn(20.));
    }
    SECTION(
        "with start time, exposure duration, and pulse duration but no PRF") {
      Waveforms::AddIsOnCalc<
          Waveforms::AddExposureDuration<Waveforms::AddStartTime<
              Waveforms::AddPulseDuration<Waveforms::Base<double>>>>>
          exp;

      exp.setExposureDuration(10);
      exp.setStartTime(2);
      exp.setPulseDuration(0.01);

      CHECK(!exp.isOn(0.0));
      CHECK(!exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(!exp.isOn(10.0));
      CHECK(!exp.isOn(11.0));
      CHECK(!exp.isOn(12.0));
      CHECK(!exp.isOn(13.0));
      CHECK(!exp.isOn(20.));
    }
    SECTION(
        "with pulse duration and PRF, but no start time or exposure duration") {
      Waveforms::AddIsOnCalc<Waveforms::AddPulseDuration<
          Waveforms::AddPulseRepititionFrequency<Waveforms::Base<double>>>>
          exp;

      exp.setPulseDuration(0.01);
      exp.setPulseRepititionFrequency(10);

      CHECK(exp.isOn(0.0));
      CHECK(exp.isOn(1.0));
      CHECK(exp.isOn(2.0));
      CHECK(exp.isOn(10.0));
      CHECK(exp.isOn(11.0));
      CHECK(exp.isOn(12.0));
      CHECK(exp.isOn(13.0));
      CHECK(exp.isOn(20.));

      CHECK(exp.isOn(0.01));
      CHECK(!exp.isOn(0.011));
      CHECK(exp.isOn(1.01));
      CHECK(!exp.isOn(1.011));
      CHECK(exp.isOn(2.01));
      CHECK(!exp.isOn(2.011));
    }
    SECTION(
        "with pulse duration, PRF, and exposure duration, but no start time") {}
    SECTION(
        "with pulse duration, PRF, start time, and exposure duration, but no "
        "start time") {}
  }
}
