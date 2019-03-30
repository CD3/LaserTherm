#include <boost/property_tree/ptree.hpp>

#include <UnitConvert.hpp>

namespace
{
namespace bpt = boost::property_tree;
}

/**
 * Converts all quantities in a property tree to values expressed
 * in the base units of a unit registry.
 */
void convertPropertyTreeUnits(bpt::ptree& config, const UnitRegistry& ureg)
{
  for (auto& it : config) {
    auto val = it.second.get_value_optional<std::string>();
    if (val) {
      try {
        auto q = ureg.makeQuantity<double>(val.value());
        it.second.put_value(q.to_base_units().value());
      } catch (...) {
      }
    }
    convertPropertyTreeUnits(it.second, ureg);
  }
}

namespace detail
{
void get_paths_imp(const bpt::ptree&                   tree,
                   std::vector<bpt::ptree::path_type>& paths,
                   bpt::ptree::path_type               current_path)
{
  if (!tree.data().empty()) {
    paths.push_back(current_path);
  }

  for (auto& it : tree) {
    auto path = current_path;
    path /= bpt::ptree::path_type(it.first, current_path.separator());
    get_paths_imp(it.second, paths, path);
  }
}
}  // namespace detail

/**
 * Return a list of all property tree paths that have data.
 */
std::vector<bpt::ptree::path_type> get_paths(const bpt::ptree& tree,
                                             char              delim = '.')
{
  std::vector<bpt::ptree::path_type> paths;
  detail::get_paths_imp(tree, paths, bpt::ptree::path_type(delim));
  return paths;
}

/**
 * Return a flattened version of the property tree.
 * All nodes are a direct child of the root/top of the tree.
 */
bpt::ptree flatten_ptree(const bpt::ptree& tree, char delim = '.')
{
  auto       paths = get_paths(tree, delim);
  bpt::ptree ftree;

  for (auto& p : paths) {
    ftree.push_back(bpt::ptree::value_type(p.dump(), tree.get_child(p).data()));
  }

  return ftree;
}

/**
 * Return an un-flattened version of the property tree. Nodes that contain the
 * delimiter in their name will be split into a node with a child node.
 *
 * @param odelim the output delimiter. the character that will be used to
 * determine the new nesting level of each node.
 * @param idelim the input delimiter. the character that is used to delimit
 * levels in the input tree.
 *
 */
bpt::ptree unflatten_ptree(const bpt::ptree& ftree, char odelim = '.',
                           char idelim = '.')
{
  auto       paths = get_paths(ftree, idelim);
  bpt::ptree tree;

  for (auto& p : paths) {
    auto key = p.dump();
    std::replace(key.begin(), key.end(), idelim, odelim);
    bpt::ptree::path_type newp(key, odelim);
    tree.put_child(newp, bpt::ptree()).put_value(ftree.get_child(p).data());
  }

  return tree;
}
