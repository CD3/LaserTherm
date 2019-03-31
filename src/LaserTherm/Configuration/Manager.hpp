#ifndef LaserTherm_Configuration_Manager_hpp
#define LaserTherm_Configuration_Manager_hpp

/** @file Manager.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/30/19
 */

#include <string>
#include <typeinfo>

#include "ptree_utils.hpp"

namespace Configuration
{
class Manager
{
 public:
  using ptree  = boost::property_tree::ptree;
  using path_t = ptree::path_type;

  ptree        configuration;
  UnitRegistry unit_registry;
  // ptree        property_units;

  template<typename T>
  T get(std::string key, std::string unit) const;
  template<typename T>
  T get(std::string key) const;

  void load(const std::string& filename);
  void load(const std::istream& in);
};
template<typename T>
T Manager::get(std::string key) const
{
  T    value;
  auto node = configuration.get_child_optional(key);
  if (!node) {
    throw std::runtime_error("A parameter named '" + key +
                             "' was requested, but was found in the tree.\nDid "
                             "you spell it correctly?");
  }

  try {
    value = node.value().get_value<T>();
  } catch (std::runtime_error& e) {
    std::stringstream error;

    error << "There was a problem converting a parameter '" << key
          << "' to the requested type.\n";
    error << "This can happen if the parameter is expected to be a numerical "
             "value, but contains non-numerical data.";
    error << "The typeid of the requested type is '" << typeid(T).name()
          << "'. The stored value for '" << key << "' is '"
          << node.value().data() << "'\n";
    error << "Do you expect the value to be convertible to the type?\n";
    error << "The error thrown during the construction is show below";
    error << "\n==========\n";
    error << e.what();
    error << "\n==========\n";

    throw std::runtime_error(error.str());
  }

  return value;
}

template<typename T>
T Manager::get(std::string key, std::string unit) const
{
  auto        raw_value = get<ptree::data_type>(key);
  Quantity<T> raw_quantity;
  T           value;
  try {
    raw_quantity = unit_registry.makeQuantity<T>(raw_value);
  } catch (std::runtime_error& e) {
    std::stringstream error;

    error << "Could not create a quantity from the requested parameter. This "
             "most likely means that the string could not be parsed as a "
             "quantity\n"
             "Do you expect this parameter to be a quantity?\n";
    error << "Parameter Name: " << key << "\n";
    error << "Parameter Value: " << raw_value << "\n";
    error << "Requested Units: " << unit << "\n";
    error << "The error thrown during the construction is show below";
    error << "\n==========\n";
    error << e.what();
    error << "\n==========\n";
    throw std::runtime_error(error.str());
  }

  try {
    value = raw_quantity.to(unit).value();
  } catch (std::runtime_error& e) {
    std::stringstream error;

    error
        << "Could not convert the requested parameter to the requested units.\n"
           "This most likely means that the parameter not have the same "
           "dimensions as the requested units.\n"
           "Did you type the parameter units correctly?\n";
    error << "Parameter Name: " << key << "\n";
    error << "Parameter Value: " << raw_value << "\n";
    error << "Requested Units: " << unit << "\n";
    error << "\n==========\n";
    error << e.what();
    error << "\n==========\n";
    throw std::runtime_error(error.str());
  }

  return value;
}
}  // namespace Configuration

#endif  // include protector
