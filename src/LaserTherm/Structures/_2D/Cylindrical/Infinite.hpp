#pragma once

/** @file Infinite.hpp
 * @brief
 * @author C.D. Clark III
 */

namespace Structures::_2D::Cylindrical
{
template<typename REAL>
class Infinite
{
 public:
  bool isInside(REAL r, REAL z) { return true; }
};

}  // namespace Structures::_2D::Cylindrical

#endif  // include protector
