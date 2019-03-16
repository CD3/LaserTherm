#include "catch.hpp"

#include<iostream>
#include<fstream>

#include <LaserTherm/Structures/_1D/Slab.hpp>
#include <LaserTherm/Structures/_1D/Collection.hpp>


TEST_CASE("Slabs")
{
  Collections::_1D::Slab< Materials::Basic<double> > slab;
  Materials::Basic<double> mat;
  mat.setConductivity( 2.0 );

  slab.setMinSurface(1.0);
  slab.setMaxSurface(2.0);

  slab.setMaterial(mat);


  Field<1,double> field(10);
  field.setCoordinateSystem( Uniform(0,5) );


  slab.set<Conductivity>( field );
  slab.set<SpecificHeat>( field );
  slab.set<AbsorptionCoefficient>( field );
  
}

TEST_CASE("Collection of Slabs")
{
}
