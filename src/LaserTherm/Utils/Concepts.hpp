#ifndef LaserTherm_Utils_Concepts_hpp
#define LaserTherm_Utils_Concepts_hpp

/** @file Concepts.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/23/19
  */

#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/member.hpp>
#include <boost/type_erasure/callable.hpp>
#include <boost/mpl/vector.hpp>

BOOST_TYPE_ERASURE_MEMBER( (has_isInside), isInside);

namespace Structures::_1D {
  template<typename T, typename Sig = bool(T)>
  using is_structure = boost::mpl::vector< has_isInside< Sig > >;
}





#endif // include protector
