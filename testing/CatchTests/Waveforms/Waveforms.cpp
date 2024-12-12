#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>

#include <LaserTherm/Waveforms/ContinuousWave.hpp>

TEST_CASE("Waveform Classes")
{
  SECTION("CW")
  {
    Waveforms::ContinuousWave<double> exp;

    exp.setExposureDuration(10);
    exp.setStartTime(2);

    CHECK(exp.getExposureDuration() == Approx(10));
    CHECK(exp.getStartTime() == Approx(2));

    CHECK(!exp.isOn(0.0));
    CHECK(!exp.isOn(1.0));
    CHECK(exp.isOn(2.0));
    CHECK(exp.isOn(10.0));
    CHECK(exp.isOn(11.0));
    CHECK(exp.isOn(12.0));
    CHECK(!exp.isOn(13.0));
    CHECK(!exp.isOn(20.));

    CHECK(exp.getWaveformValue(0.0) == Approx(0));
    CHECK(exp.getWaveformValue(1.0) == Approx(0));
    CHECK(exp.getWaveformValue(2.0) == Approx(1));
    CHECK(exp.getWaveformValue(10.0) == Approx(1));
    CHECK(exp.getWaveformValue(11.0) == Approx(1));
    CHECK(exp.getWaveformValue(12.0) == Approx(1));
    CHECK(exp.getWaveformValue(13.0) == Approx(0));
    CHECK(exp.getWaveformValue(20.) == Approx(0));
  }
}
