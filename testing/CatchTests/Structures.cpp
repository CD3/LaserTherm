#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <Field.hpp>
#include <LaserTherm/Materials/Basic.hpp>
#include <LaserTherm/Structures/_1D/Slab.hpp>
//#include <LaserTherm/Structures/_1D/Collection.hpp>

TEST_CASE("Slabs")
{
  Structures::_1D::Slab<Materials::Basic<double>, double> slab;
  Materials::Basic<double>                                mat;
  mat.setConductivity(2.0);

  SECTION("1 surface")
  {
    SECTION("min surface")
    {
      slab.setMinSurfacePosition(1.0);

      CHECK(!slab.isInside(0));
      CHECK(!slab.isInside(0.5));
      CHECK( slab.isInside(1));
      CHECK( slab.isInside(1.5));
      CHECK( slab.isInside(2.0));
      CHECK( slab.isInside(2.5));

      CHECK(!slab.isInside(0.99));
      CHECK( slab.isInside(1.01));
      CHECK( slab.isInside(1.99));
      CHECK( slab.isInside(2.01));
    }
    SECTION("max surface")
    {
      slab.setMaxSurfacePosition(1.0);

      CHECK( slab.isInside(0));
      CHECK( slab.isInside(0.5));
      CHECK( slab.isInside(1));
      CHECK(!slab.isInside(1.5));
      CHECK(!slab.isInside(2.0));
      CHECK(!slab.isInside(2.5));

      CHECK( slab.isInside(0.99));
      CHECK(!slab.isInside(1.01));
      CHECK(!slab.isInside(1.99));
      CHECK(!slab.isInside(2.01));
    }
  }

  SECTION("2 surface")
  {
    slab.setMinSurfacePosition(1.0);
    slab.setMaxSurfacePosition(2.0);

    CHECK(!slab.isInside(0));
    CHECK(!slab.isInside(0.5));
    CHECK(slab.isInside(1));
    CHECK(slab.isInside(1.5));
    CHECK(slab.isInside(2.0));
    CHECK(!slab.isInside(2.5));

    CHECK(!slab.isInside(0.99));
    CHECK(slab.isInside(1.01));
    CHECK(slab.isInside(1.99));
    CHECK(!slab.isInside(2.01));
  }

  // slab.setMaterial(mat);

  Field<double, 1> field(10);
  field.setCoordinateSystem(Uniform(0., 5.));

  // set<Conductivity>( field );
  // set<SpecificHeat>( field );
  // set<AbsorptionCoefficient>( field );

  // slab.set( field, Materials::Basic<double>::getConductivity );
}

TEST_CASE("Collection of Slabs") {}
