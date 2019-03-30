#include "catch.hpp"

/**
 * This file is used for developement. As new classes are created, small tests
 * are written here so that we can try to compile and use them.
 */


#include<UnitConvert.hpp>

TEST_CASE("General", "[devel]")
{
  UnitRegistry ureg;

  ureg.addBaseUnit<Dimension::Name::Length>("cm");
  ureg.addBaseUnit<Dimension::Name::Mass>("g");
  ureg.addBaseUnit<Dimension::Name::Time>("s");
  ureg.addBaseUnit<Dimension::Name::Temperature>("K");
  ureg.addBaseUnit<Dimension::Name::Amount>("mol");
  
}
