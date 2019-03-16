#include "catch.hpp"

#include <iostream>
#include <fstream>

#include <LaserTherm/Materials/Basic.hpp>

TEST_CASE("Basic Material Tests")
{

  Materials::Basic<double> mat1, mat2;

  mat1.setDensity(1);
  mat1.setSpecificHeat(2);
  mat1.setConductivity(3);

  mat2.setDensity(10);
  mat2.setSpecificHeat(20);
  mat2.setConductivity(30);
  mat2.setAbsorptionCoefficient(100);
  mat2.setScatteringCoefficient(200);


  CHECK( (bool)mat1.getDensity() );
  CHECK( (bool)mat1.getSpecificHeat() );
  CHECK( (bool)mat1.getConductivity() );
  CHECK(!(bool)mat1.getAbsorptionCoefficient() );
  CHECK(!(bool)mat1.getScatteringCoefficient() );

  CHECK( mat1.getDensity() == Approx(1 ) );
  CHECK( mat2.getDensity() == Approx(10) );


  CHECK( (bool)mat2.getDensity() );
  CHECK( (bool)mat2.getSpecificHeat() );
  CHECK( (bool)mat2.getConductivity() );
  CHECK( (bool)mat2.getAbsorptionCoefficient() );
  CHECK( (bool)mat2.getScatteringCoefficient() );


}
