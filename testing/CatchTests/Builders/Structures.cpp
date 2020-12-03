#include "catch.hpp"

#include <fstream>
#include <iostream>
#include <optional>

#include <boost/property_tree/ptree.hpp>

#include <libField/Field.hpp>
#include <LaserTherm/MaterialStructure.hpp>
#include <LaserTherm/Materials/Basic.hpp>
#include <LaserTherm/Structures/_1D/AnyStructure.hpp>
#include <LaserTherm/Configuration/Builders.hpp>

TEST_CASE("Material Map Builder")
{
  boost::property_tree::ptree config;

  config.put("simulation.dimensions", 1);

  config.put("layers.0.material", "water");

  config.put("materials.water.thermal.density", 1 );
  config.put("materials.water.thermal.conductivity", 0.00628 );
  config.put("materials.water.thermal.specific_heat", 4.1868 );
  config.put("materials.water.optical.absorption_coefficient", 10 );
  config.put("materials.water.optical.scattering_coefficient", 1 );
  config.put("materials.water.optical.anisotropy", 0.9 );
  config.put("materials.water.thermal.bc.convection.transfer_rate", 1e-3 );

  config.put("materials.steel.thermal.density", 8.05 );
  config.put("materials.steel.thermal.conductivity", 0.45 );
  config.put("materials.steel.thermal.specific_heat", 0.49 );

  config.put("layers.1.description", "absorber" );
  config.put("layers.1.optical.absorption_coefficient", 1 );
  config.put("layers.1.position", 0 );
  config.put("layers.1.thickness", 10e-4 );


  std::map< std::string, Materials::Basic<double> > materials;
  Builders::build(materials, config.get_child("materials"));

  CHECK( materials.size() == 2 );

  REQUIRE( materials.find("water") != materials.end() );
  REQUIRE( materials.find("steel") != materials.end() );

  CHECK( materials["water"].getSpecificHeatCapacity() );
  CHECK( materials["water"].getDensity() );
  CHECK( materials["water"].getThermalConductivity() );
  CHECK( materials["water"].getAbsorptionCoefficient() );
  CHECK( materials["water"].getScatteringCoefficient() );

  CHECK( materials["steel"].getSpecificHeatCapacity() );
  CHECK( materials["steel"].getDensity() );
  CHECK( materials["steel"].getThermalConductivity() );
  CHECK(!materials["steel"].getAbsorptionCoefficient() );
  CHECK(!materials["steel"].getScatteringCoefficient() );

  CHECK( materials["water"].getSpecificHeatCapacity().value() == Approx(4.1868) );
  CHECK( materials["water"].getAbsorptionCoefficient().value() == Approx(10) );


}

TEST_CASE("1D Structure Set Builder")
{
  boost::property_tree::ptree config;

  config.put("simulation.dimensions", 1);


  config.put("materials.water.thermal.density", 1 );
  config.put("materials.water.thermal.conductivity", 0.00628 );
  config.put("materials.water.thermal.specific_heat", 4.1868 );
  config.put("materials.water.thermal.bc.convection.transfer_rate", 1e-3 );

  config.put("layers.0.material", "water");

  config.put("layers.1.description", "absorber");
  config.put("layers.1.optical.absorption_coefficient", 10 );
  config.put("layers.1.position", 0.1 );
  config.put("layers.1.thickness", 0.2 );

  config.put("layers.2.description", "absorber");
  config.put("layers.2.optical.absorption_coefficient", 100 );
  config.put("layers.2.thickness", 0.3 );

  std::vector<
      MaterialStructure<Materials::Basic<double>, Structures::_1D::AnyStructure<double> > >
      structures;

  Builders::build(structures, config);

  Field<double,1> f(11);
  f.setCoordinateSystem(Uniform(0,1));

  // | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |  position
  // |   0 |   1 |   1 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |   0 |  material

  f.set(0);
  for (auto& s : structures) {
    f.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getThermalConductivity())
        return s.material.getThermalConductivity().value();
      return std::nullopt;
    });
  }

  CHECK( f( 0) == Approx(0.00628) );
  CHECK( f( 1) == Approx(0.00628) );
  CHECK( f( 2) == Approx(0.00628) );
  CHECK( f( 3) == Approx(0.00628) );
  CHECK( f( 4) == Approx(0.00628) );
  CHECK( f( 5) == Approx(0.00628) );
  CHECK( f( 6) == Approx(0.00628) );
  CHECK( f( 7) == Approx(0.00628) );
  CHECK( f( 8) == Approx(0.00628) );
  CHECK( f( 9) == Approx(0.00628) );
  CHECK( f(10) == Approx(0.00628) );

  // | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |  position
  // |   0 |  10 |  10 | 100 | 100 | 100 |   0 |   0 |   0 |   0 |   0 |  material
  f.set(0);
  for (auto& s : structures) {
    f.set_f([&s](auto x) -> std::optional<double> {
        auto mua = s.material.getAbsorptionCoefficient();
      if (s.structure.isInside(x[0]) && s.material.getAbsorptionCoefficient())
        return s.material.getAbsorptionCoefficient().value();
      return std::nullopt;
    });
  }


  CHECK( f( 0) == Approx(0  ) );
  CHECK( f( 1) == Approx(0  ) );
  CHECK( f( 2) == Approx(10 ) );
  CHECK( f( 3) == Approx(10 ) );
  CHECK( f( 4) == Approx(100) );
  CHECK( f( 5) == Approx(100) );
  CHECK( f( 6) == Approx(100) );
  CHECK( f( 7) == Approx(0  ) );
  CHECK( f( 8) == Approx(0  ) );
  CHECK( f( 9) == Approx(0  ) );
  CHECK( f(10) == Approx(0  ) );

}
