#include "catch.hpp"

#include <iostream>

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <LaserTherm/Configuration/ptree_utils.hpp>

TEST_CASE("Property Tree Unit Conversion")
{
  UnitConvert::UnitRegistry ureg;

  ureg.addBaseUnit<UnitConvert::Dimension::Name::Length>("cm");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Mass>("g");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Time>("s");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Temperature>("K");
  ureg.addBaseUnit<UnitConvert::Dimension::Name::Amount>("mol");

  ureg.addUnit("m = 100 cm");
  ureg.addUnit("L = 1000 cm^3");
  ureg.addUnit("J = kg m^2 / s^2");
  ureg.addUnit("cal = 4.18 J");

  boost::property_tree::ptree config;

  config.put("simulation.grid.x.min", "-10 cm");
  config.put("simulation.grid.x.max", "1 m");

  config.put("layers.0.thermal.density", "1 g/mL");
  config.put("layers.0.thermal.specific_heat", "1 cal/g/K");
  config.put("layers.0", "water");

  CHECK(config.get<std::string>("layers.0") == "water");
  CHECK(config.get<std::string>("layers.0.thermal.density") == "1 g/mL");

  convertPropertyTreeUnits(config, ureg);

  boost::property_tree::write_xml(std::cout, config);
}

TEST_CASE("Property Tree Utils")
{
  SECTION("Get all paths in tree")
  {
    using path_type = boost::property_tree::ptree::path_type;

    boost::property_tree::ptree tree;
    tree.put("l1.l1.l1", "111");
    tree.put("l1.l1.l2", "112");
    tree.put("l1.l1", "11");
    tree.put("l2.l1.l1", "211");
    tree.put("", "0");
    tree.put(path_type("l1.0|l1.0", '|'), "3010");

    auto paths = get_paths(tree, '|');
    CHECK(paths.size() == 6);
    CHECK(paths[0].dump() == "");
    CHECK(paths[1].dump() == "l1|l1");
    CHECK(paths[5].dump() == "l1.0|l1.0");
  }

  SECTION("Flatten tree")
  {
    boost::property_tree::ptree tree;
    tree.put("a1.b1.c1", "111");
    tree.put("a1.b1.c2", "112");
    tree.put("a2.b1", "21");
    tree.put("a3.b1", "31");
    tree.put("a3.b1.c1", "311");
    tree.put(boost::property_tree::ptree::path_type("a4|b1.1", '|'), "41.1");

    SECTION("Unique Delimiter Char")
    {
      using path_t = boost::property_tree::ptree::path_type;
      auto ftree   = flatten_ptree(tree, '|');
      CHECK(tree.size() == 4);
      CHECK(ftree.size() == 6);
      CHECK(tree.get_child("a1.b1.c1").data() ==
            ftree.get_child("a1|b1|c1").data());
      CHECK(tree.get_child("a1.b1.c2").data() ==
            ftree.get_child("a1|b1|c2").data());
      CHECK(tree.get_child("a2.b1").data() == ftree.get_child("a2|b1").data());
      CHECK(tree.get_child("a3.b1").data() == ftree.get_child("a3|b1").data());
      CHECK(tree.get_child("a3.b1.c1").data() ==
            ftree.get_child("a3|b1|c1").data());
      CHECK(tree.get_child(path_t("a4/b1.1", '/')).data() ==
            ftree.get_child(path_t("a4|b1.1", '/')).data());
    }
  }

  SECTION("Un-flatten tree")
  {
    boost::property_tree::ptree ftree;
    ftree.put("a1|b1|c1", "111");
    ftree.put("a1|b1|c2", "112");
    ftree.put("a2|b1", "21");
    ftree.put("a3|b1", "31");
    ftree.put("a3|b1.c1", "311");  // a two-level node

    using path_t = boost::property_tree::ptree::path_type;
    auto tree    = unflatten_ptree(ftree, '|');
    CHECK(ftree.size() == 4);
    CHECK(tree.size() == 3);
    CHECK(tree.get_child("a1.b1.c1").data() == ftree.get_child("a1|b1|c1").data());
    CHECK(tree.get_child("a3.b1.c1").data() == ftree.get_child("a3|b1.c1").data());
  }

  SECTION("Custom Integer Translator")
  {
    boost::property_tree::ptree tree;
    tree.put("x", "1");
    tree.put("y", "1.0");

    CHECK( tree.get<double>("x") == Approx(1) );
    CHECK( tree.get<double>("y") == Approx(1) );
    CHECK( tree.get<int>("x") == 1 );
    CHECK( tree.get<int>("y") == 1 );


  }
}
