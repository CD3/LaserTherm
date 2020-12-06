#include "catch.hpp"

#include <fstream>
#include <iostream>

#include <boost/property_tree/info_parser.hpp>

#include <LaserTherm/Configuration/Manager.hpp>

TEST_CASE("Configuration Manager", "[configuration]")
{
  std::string config_text = R"(
  simulation.dimensions = 1
  simulation.grid.x.min = 0
  simulation.grid.x.max = 5 cm
  layers.0.material = water
  layers.0.position = 0
  layers.0.thickness = 5 mm
  sensors.0.type = temperature
  simulation.grid.z.max = 10 inch
  )";

  Configuration::Manager config;

  SECTION("Loading")
  {
    SECTION("from INI file")
    {
      std::stringstream in(config_text);

      std::ofstream out("config-test.ini");
      out << config_text;
      out.close();

      config.load("config-test.ini");

      CHECK(config.configuration.size() == 3);
    }
  }

  SECTION("Element Access")
  {
    std::ofstream out("config-test.ini");
    out << config_text;
    out.close();

    config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Length>(
        "cm");
    config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Mass>("g");
    config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Time>("s");
    config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Temperature>(
        "K");
    config.unit_registry.addBaseUnit<UnitConvert::Dimension::Name::Amount>(
        "mol");
    config.unit_registry
        .addBaseUnit<UnitConvert::Dimension::Name::ElectricalCurrent>("A");
    config.unit_registry
        .addBaseUnit<UnitConvert::Dimension::Name::LuminousIntensity>("cd");

    config.unit_registry.addUnit("m = 100 cm");
    config.unit_registry.addUnit("L = 1000 cm^3");
    config.unit_registry.addUnit("J = kg m^2 / s^2");
    config.unit_registry.addUnit("W = J/s");
    config.unit_registry.addUnit("cal = 4.184 J");
    config.unit_registry.addUnit("degC = K - 273.15");

    config.load("config-test.ini");

    CHECK(config.configuration.get<std::string>("simulation.dimensions") ==
          "1");

    SECTION("Access Errors")
    {
      CHECK_THROWS_WITH(config.get<std::string>("missing"),
                        Catch::StartsWith("A parameter named 'missing'"));
      CHECK_THROWS_WITH(config.get<double>("simulation.grid.x.max"),
                        Catch::StartsWith("There was a problem converting") &&
                            Catch::Contains("simulation.grid.x.max"));
      CHECK_THROWS_WITH(config.get_quantity<double>("layers.0.material", "cm"),
                        Catch::Contains("not create a quantity"));
      CHECK_THROWS_WITH(
          config.get_quantity<double>("simulation.grid.x.max", "s"),
          Catch::Contains("not convert the requested parameter") &&
              Catch::Contains("Dimension Error"));
    }

    CHECK(config.get<std::string>("simulation.dimensions") == "1");
    CHECK(config.get<int>("simulation.dimensions") == 1);
    CHECK(config.get_quantity<double>("simulation.grid.x.max", "cm") ==
          Approx(5));
    CHECK(config.get_quantity<double>("simulation.grid.x.max", "m") ==
          Approx(0.05));

    SECTION("With Root Node")
    {
      config.root = Configuration::Manager::path_t("layers.0");

      CHECK_THROWS_WITH(config.get<std::string>("missing"),
                        Catch::StartsWith("A parameter named 'missing'"));
      CHECK(config.get<std::string>("material") == "water");
      CHECK(config.get_quantity<double>("thickness", "cm") == Approx(0.5));
    }

    SECTION("Query Element Existence")
    {
      CHECK(config.has("simulation"));
      CHECK(config.has("layers.0"));
      CHECK(config.has("layers.0.material"));
      CHECK(!config.has("missing"));
      CHECK(!config.has("layers.0.missing"));

      config.root = Configuration::Manager::path_t("layers.0");
      CHECK(!config.has("simulation"));
      CHECK(config.has("material"));
    }
  }
}
