
#include "./Manager.hpp"

#include <fstream>
#include <sstream>

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace Configuration
{
void Manager::load(const std::string& filename)
{
  std::stringstream error;
  error << "Could not parse configuration file: " << filename << "\n";
  bool fail = true;

  if (fail) {  // INI
    std::ifstream in(filename);
    try {
      boost::property_tree::read_ini(in, this->configuration);
      fail = false;
    } catch (boost::property_tree::ini_parser_error& e) {
      error
          << "\n\nWhen trying to parse file as INI, parser threw this error:\n";
      error << "\n==========\n";
      error << e.what();
      error << "\n==========\n";
    }
  }

  if (fail) {  // JSON
    std::ifstream in(filename);
    try {
      boost::property_tree::read_json(in, this->configuration);
      fail = false;
    } catch (boost::property_tree::json_parser_error& e) {
      error << "\n\nWhen trying to parse file as JSON, parser threw this "
               "error:\n";
      error << "\n==========\n";
      error << e.what();
      error << "\n==========\n";
    }
  }

  if (fail) {  // XML
    std::ifstream in(filename);
    try {
      boost::property_tree::read_xml(in, this->configuration);
      fail = false;
    } catch (boost::property_tree::xml_parser_error& e) {
      error
          << "\n\nWhen trying to parse file as XML, parser threw this error:\n";
      error << "\n==========\n";
      error << e.what();
      error << "\n==========\n";
    }
  }

  if (fail) {  // INFO
    std::ifstream in(filename);
    try {
      boost::property_tree::read_info(in, this->configuration);
      fail = false;
    } catch (boost::property_tree::info_parser_error& e) {
      error << "\n\nWhen trying to parse file as INFO, parser threw this "
               "error:\n";
      error << "\n==========\n";
      error << e.what();
      error << "\n==========\n";
    }
  }

  if (fail) throw std::runtime_error(error.str());

  this->configuration = unflatten_ptree(this->configuration);

  return;
}

bool Manager::has(const path_t& path) const
{
  return static_cast<bool>(configuration.get_child_optional(addRoot(path)));
}

Manager::path_t Manager::addRoot(const path_t& path) const
{
  return root / path;
}
}  // namespace Configuration
