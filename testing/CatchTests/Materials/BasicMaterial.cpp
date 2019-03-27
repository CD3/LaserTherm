#include "catch.hpp"

#include <iostream>
#include <fstream>

#include <LaserTherm/Materials/Basic.hpp>

TEST_CASE("Basic Material Tests")
{

  Materials::Basic<double> mat1, mat2;

  mat1.setDensity(1);
  mat1.setSpecificHeatCapacity(2);
  mat1.setThermalConductivity(3);

  mat2.setDensity(10);
  mat2.setSpecificHeatCapacity(20);
  mat2.setThermalConductivity(30);
  mat2.setAbsorptionCoefficient(100);
  mat2.setScatteringCoefficient(200);


  CHECK( (bool)mat1.getDensity() );
  CHECK( (bool)mat1.getSpecificHeatCapacity() );
  CHECK( (bool)mat1.getThermalConductivity() );
  CHECK(!(bool)mat1.getAbsorptionCoefficient() );
  CHECK(!(bool)mat1.getScatteringCoefficient() );

  CHECK( mat1.getDensity() == Approx(1 ) );
  CHECK( mat2.getDensity() == Approx(10) );


  CHECK( (bool)mat2.getDensity() );
  CHECK( (bool)mat2.getSpecificHeatCapacity() );
  CHECK( (bool)mat2.getThermalConductivity() );
  CHECK( (bool)mat2.getAbsorptionCoefficient() );
  CHECK( (bool)mat2.getScatteringCoefficient() );


}
