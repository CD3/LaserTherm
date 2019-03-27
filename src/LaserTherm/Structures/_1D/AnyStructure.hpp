#ifndef LaserTherm_Structures__1D_AnyStructure_hpp
#define LaserTherm_Structures__1D_AnyStructure_hpp

/** @file AnyStructure.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/18/19
  */

#include "../../Utils/Concepts.hpp"

namespace Structures::_1D
{



template<typename T, typename Sig = bool(T)>
using AnyStructure = boost::type_erasure::any< boost::mpl::vector<
                                                       boost::type_erasure::copy_constructible<>,
                                                       is_structure<T,Sig>,
                                                       boost::type_erasure::relaxed >, boost::type_erasure::_self >;
}


#endif // include protector
