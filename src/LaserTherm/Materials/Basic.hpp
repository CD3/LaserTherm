#ifndef LaserTherm_Materials_BasicMaterial_hpp
#define LaserTherm_Materials_BasicMaterial_hpp

/** @file BasicMaterial.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/02/19
 */
#include "./Mixins.hpp"

namespace Materials
{
template<typename REAL>
using Basic = AddOptionalSpecificHeatCapacity<
    AddOptionalDensity<AddOptionalThermalConductivity<
        AddOptionalAbsorptionCoefficient<AddOptionalScatteringCoefficient<
            AddOptionalAnisotropy<Base<REAL> > > > > > >;

}

#endif  // include protector
