#ifndef LaserTherm_Structures__1D_Infinite_hpp
#define LaserTherm_Structures__1D_Infinite_hpp

/** @file Infinite.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/23/19
 */

namespace Structures::_1D
{
template<typename REAL>
class Infinite
{
 public:
  bool isInside(REAL x) { return true; }
};

}  // namespace Structures::_1D

#endif  // include protector
