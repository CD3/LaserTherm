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

MAKE_ADD_OPTIONAL_MEMBER_MIXIN(SpecificHeatCapacity);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(Density);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(ThermalConductivity);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(VolumetricSpecificHeatCapacity);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(AbsorptionCoefficient);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(ScatteringCoefficient);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(Anisotropy);
MAKE_ADD_OPTIONAL_MEMBER_MIXIN(RefractiveIndex);

}  // namespace Materials

#endif  // include protector
