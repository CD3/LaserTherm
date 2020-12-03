#ifndef LaserTherm_Utils_TypeTraits_hpp
#define LaserTherm_Utils_TypeTraits_hpp

/** @file TypeTraits.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/10/19
  */

#include <libField/Field.hpp>

template<typename T>
struct get_real_type {};

template<typename REAL>
struct get_real_type<Field<REAL,1>> { using type = REAL; };



#endif // include protector
