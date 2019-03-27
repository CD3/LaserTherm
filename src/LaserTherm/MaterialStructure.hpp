#ifndef LaserTherm_MaterialStructure_hpp
#define LaserTherm_MaterialStructure_hpp

/** @file MaterialStructure.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/18/19
 */


template<class MaterialType, class StructureType>
class MaterialStructure
{
 public:
  MaterialType  material;
  StructureType structure;

  MaterialStructure( MaterialType m, StructureType s) : structure(s), material(m)
  {
  }
  MaterialStructure() = default;

}; // namespace Structures::_1D

#endif  // include protector
