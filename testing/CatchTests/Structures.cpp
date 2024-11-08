#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
using namespace Catch;

#include <fstream>
#include <iostream>

#include <boost/optional.hpp>

#include <LaserTherm/MaterialStructure.hpp>
#include <LaserTherm/Materials/Basic.hpp>
#include <LaserTherm/Structures/_1D/AnyStructure.hpp>
#include <LaserTherm/Structures/_1D/Cartesian/Slab.hpp>
#include <LaserTherm/Structures/_1D/Infinite.hpp>
#include <libField/Field.hpp>

TEST_CASE("Structures") {
  SECTION("Slab") {
    Structures::_1D::Cartesian::Slab<double> slab;

    SECTION("1 surface") {
      SECTION("min surface") {
        slab.setMinSurfacePosition(1.0);

        CHECK(!slab.isInside(0));
        CHECK(!slab.isInside(0.5));
        CHECK(slab.isInside(1));
        CHECK(slab.isInside(1.5));
        CHECK(slab.isInside(2.0));
        CHECK(slab.isInside(2.5));

        CHECK(!slab.isInside(0.99));
        CHECK(slab.isInside(1.01));
        CHECK(slab.isInside(1.99));
        CHECK(slab.isInside(2.01));
      }
      SECTION("max surface") {
        slab.setMaxSurfacePosition(1.0);

        CHECK(slab.isInside(0));
        CHECK(slab.isInside(0.5));
        CHECK(slab.isInside(1));
        CHECK(!slab.isInside(1.5));
        CHECK(!slab.isInside(2.0));
        CHECK(!slab.isInside(2.5));

        CHECK(slab.isInside(0.99));
        CHECK(!slab.isInside(1.01));
        CHECK(!slab.isInside(1.99));
        CHECK(!slab.isInside(2.01));
      }
    }

    SECTION("2 surface") {
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

    Field<double, 1> field(10);
    field.setCoordinateSystem(Uniform(0., 5.));
  }
  SECTION("Infinite") {
    Structures::_1D::Infinite<double> structure;

    CHECK(structure.isInside(-100));
    CHECK(structure.isInside(-1));
    CHECK(structure.isInside(0));
    CHECK(structure.isInside(1));
    CHECK(structure.isInside(100));
  }
  SECTION("AnyStructure") {
    Structures::_1D::AnyStructure<double> structure;
    Structures::_1D::Cartesian::Slab<double> slab;
    slab.setMinSurfacePosition(0);
    slab.setThickness(1.5);
    structure = slab;

    CHECK(!structure.isInside(-1));
    CHECK(!structure.isInside(-0.1));
    CHECK(structure.isInside(0));
    CHECK(structure.isInside(1.0));
    CHECK(structure.isInside(1.5));
    CHECK(!structure.isInside(1.6));

    structure = Structures::_1D::Infinite<double>();

    CHECK(structure.isInside(-1));
    CHECK(structure.isInside(-0.1));
    CHECK(structure.isInside(0));
    CHECK(structure.isInside(1.0));
    CHECK(structure.isInside(1.5));
    CHECK(structure.isInside(1.6));
  }
}

TEST_CASE("Material Structure") {
  MaterialStructure<Materials::Basic<double>,
                    Structures::_1D::Cartesian::Slab<double>>
      obj;

  obj.structure.setMinSurfacePosition(1.0);
  obj.structure.setMaxSurfacePosition(2.0);
  obj.material.setThermalConductivity(2);
  obj.material.setSpecificHeatCapacity(3);
  obj.material.setDensity(4);

  Field<double, 1> f(101);
  f.setCoordinateSystem(Uniform(0, 10));

  f.set(0);

  f.set_f([&obj](auto x) -> std::optional<double> {
    if (obj.structure.isInside(x[0]))
      return obj.material.getThermalConductivity().value();
    return std::nullopt;
  });

  CHECK(f(9) == Approx(0));
  CHECK(f(10) == Approx(2));
  CHECK(f(11) == Approx(2));
  CHECK(f(20) == Approx(2));
  CHECK(f(21) == Approx(0));

  f.set_f([&obj](auto x) -> std::optional<double> {
    if (obj.structure.isInside(x[0]))
      return obj.material.getDensity().value() *
             obj.material.getSpecificHeatCapacity().value();
    return std::nullopt;
  });

  CHECK(f(9) == Approx(0));
  CHECK(f(10) == Approx(12));
  CHECK(f(11) == Approx(12));
  CHECK(f(20) == Approx(12));
  CHECK(f(21) == Approx(0));
}
TEST_CASE("Collection of Material Structures") {
  using MaterialType = Materials::Basic<double>;
  using StructureType = Structures::_1D::AnyStructure<double>;
  using MaterialStructureType = MaterialStructure<MaterialType, StructureType>;

  std::vector<MaterialStructureType> structures;

  MaterialType mat;
  mat.setThermalConductivity(2);
  mat.setDensity(1);

  structures.push_back(
      MaterialStructureType(mat, Structures::_1D::Infinite<double>()));

  mat = MaterialType();
  Structures::_1D::Cartesian::Slab<double> slab;
  for (int i = 0; i < 3; ++i) {
    MaterialType mat;
    mat.setDensity(1.5 * (i + 1));
    if (i == 0)
      slab.setMinSurfacePosition(0.0);
    else
      slab.setMinSurfacePosition(slab.getMaxSurfacePosition());

    slab.setThickness(0.5 * (i + 1));

    structures.push_back(MaterialStructureType(mat, slab));
  }

  Field<double, 1> rho(101);
  rho.setCoordinateSystem(Uniform(-1, 9));
  Field<double, 1> k(rho.getCoordinateSystemPtr());

  for (auto &s : structures) {
    rho.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getDensity())
        return s.material.getDensity().value();
      return std::nullopt;
    });
  }

  for (auto &s : structures) {
    k.set_f([&s](auto x) -> std::optional<double> {
      if (s.structure.isInside(x[0]) && s.material.getThermalConductivity())
        return s.material.getThermalConductivity().value();
      return std::nullopt;
    });
  }

  CHECK(rho(0) == Approx(1));
  CHECK(rho(9) == Approx(1));

  CHECK(rho(10) == Approx(1.5));
  CHECK(rho(14) == Approx(1.5));

  CHECK(rho(15) == Approx(3));
  CHECK(rho(20) == Approx(3));
  CHECK(rho(24) == Approx(3));

  CHECK(rho(25) == Approx(4.5));
  CHECK(rho(35) == Approx(4.5));
  CHECK(rho(39) == Approx(4.5));
  CHECK(rho(40) == Approx(4.5));

  CHECK(rho(41) == Approx(1));
  CHECK(rho(50) == Approx(1));

  CHECK(rho(100) == Approx(1));
}
