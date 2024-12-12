#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>

#include <LaserTherm/HeatSources/_1D/Cartesian/BeersLaw.hpp>

TEST_CASE("Beers Law Source")
{
  HeatSources::_1D::Cartesian::BeersLaw<double> HeatSource(11);
  HeatSource.A.setCoordinateSystem(Uniform(-2, 2));
  HeatSource.E0 = 10;

  SECTION("Single Absorption Coefficient")
  {
    HeatSource.mu_a.set(2);
    HeatSource.compute();

    CHECK(HeatSource.A(0) == Approx(2 * 10));
    CHECK(HeatSource.A(10) == Approx(2 * 10 * exp(-2 * 4)));
  }

  SECTION("Two Absorption Coefficient")
  {
    HeatSource.mu_a.set_f([](auto i, auto cs) {
      if (i[0] <= 5) return 2;
      return 4;
    });
    HeatSource.compute();
    std::ofstream out("A.txt");
    out << HeatSource.A;

    CHECK(HeatSource.A(0) == Approx(2 * 10));
    CHECK(HeatSource.A(10) == Approx(4 * 10 * exp(-2 * 2 - 4 * 2)));
  }
}
