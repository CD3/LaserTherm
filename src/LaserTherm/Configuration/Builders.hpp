#ifndef LaserTherm_Configuration_Builders_hpp
#define LaserTherm_Configuration_Builders_hpp

/** @file Builders.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/24/19
 */

#include <iostream>
#include <map>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "../MaterialStructure.hpp"
#include "../Materials/Basic.hpp"
#include "../Structures/_1D/AnyStructure.hpp"
#include "../Structures/_1D/Infinite.hpp"
#include "../Structures/_1D/Slab.hpp"

namespace Builders
{
template<class REAL>
void buildMaterial(Materials::Basic<REAL>&            material,
                   const boost::property_tree::ptree& config)
{
  boost::optional<REAL> prop;

  /* see ../Materials/Basic.hpp for list of properties */
  prop = config.get_optional<REAL>("thermal.specific_heat");
  if (prop) material.setSpecificHeatCapacity(prop.value());

  prop = config.get_optional<REAL>("thermal.density");
  if (prop) material.setDensity(prop.value());

  prop = config.get_optional<REAL>("thermal.conductivity");
  if (prop) material.setThermalConductivity(prop.value());

  prop = config.get_optional<REAL>("optical.absorption_coefficient");
  if (prop) material.setAbsorptionCoefficient(prop.value());

  prop = config.get_optional<REAL>("optical.scattering_coefficient");
  if (prop) material.setScatteringCoefficient(prop.value());
}

template<class REAL>
void buildMaterialMap(std::map<std::string, Materials::Basic<REAL>>& materials,
                      const boost::property_tree::ptree&             config)
{
  using MaterialType = Materials::Basic<REAL>;

  for (auto it : config) {
    MaterialType mat;
    buildMaterial(mat, it.second);
    materials[it.first] = mat;
  }
}

/**
 * Builds an vector of 1D material layers from a configuration defined in a
 * boost property_tree.
 */
template<class REAL>
void buildStructures(
    std::vector<MaterialStructure<Materials::Basic<REAL>,
                                  Structures::_1D::AnyStructure<REAL>>>&
                                       structures,
    const boost::property_tree::ptree& config)
{
  using MaterialType          = Materials::Basic<REAL>;
  using StructureType         = Structures::_1D::AnyStructure<REAL>;
  using MaterialStructureType = MaterialStructure<MaterialType, StructureType>;
  using InfStucture           = Structures::_1D::Infinite<REAL>;
  using LayerStucture         = Structures::_1D::Slab<REAL>;

  std::map<std::string, MaterialType> materials;
  auto materials_config = config.get_child_optional("materials");
  if (materials_config) {
    buildMaterialMap(materials, materials_config.value());
  }

  auto layers_config = config.get_child_optional("layers");
  if (layers_config) {
    boost::optional<REAL> current_position;
    for (auto it : layers_config.value()) {
      MaterialType mat;

      if (it.second.find("material") != it.second.not_found())
        mat = materials[it.second.get<std::string>("material")];
      else
        buildMaterial(mat,it.second);

      auto position  = it.second.get_optional<REAL>("position");
      auto thickness = it.second.get_optional<REAL>("thickness");

      if (!position && !thickness && !current_position) {
        structures.push_back(
            MaterialStructureType(mat, StructureType(InfStucture())));
        continue;
      }

      LayerStucture layer;
      if (!position && current_position) {
        layer.setMinSurfacePosition(current_position.value());
      }

      if (position && !current_position) {
        layer.setMinSurfacePosition(position.value());
      }

      if (thickness) {
        layer.setThickness(thickness.value());
        current_position = layer.getMaxSurfacePosition();
      }

      structures.push_back(MaterialStructureType(mat, layer));
    }
  }
}

}  // namespace Builders

#endif  // include protector
