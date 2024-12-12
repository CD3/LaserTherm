#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

/*
 * This file is used for developement. As new classes are created, small tests
 * are written here so that we can try to compile and use them.
 */

#include <UnitConvert.hpp>

TEST_CASE("General", "[devel]")
{
  UnitConvert::UnitRegistry ureg;

  ureg.addBaseUnit<UnitConvert::Dimension::Name::Length>("cm");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Mass>("g");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Time>("s");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Temperature>("K");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Amount>("mol");
}
