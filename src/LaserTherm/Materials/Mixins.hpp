#ifndef LaserTherm_Materials_Mixins_hpp
#define LaserTherm_Materials_Mixins_hpp

/** @file Mixins.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/02/19
  */

#include "../Utils/Mixins.hpp"

namespace Materials
{

  template<typename REAL>
  using Base = MixinBase<REAL>;

  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( SpecificHeat );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( Density );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( Conductivity );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( VolumetricSpecificHeat );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( AbsorptionCoefficient );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( ScatteringCoefficient );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( Anisotropy );
  MAKE_ADD_OPTIONAL_MEMBER_MIXIN( RefractiveIndex );

}




#endif // include protector
