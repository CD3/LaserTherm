#ifndef LaserTherm_Waveforms_RegularPulseTrain_hpp
#define LaserTherm_Waveforms_RegularPulseTrain_hpp

/** @file RegularPulseTrain.cpp
 * @brief
 * @author C.D. Clark III
 * @date 02/28/19
 */

namespace Waveforms
{
template<typename REAL>
using RegularPulseTrain =
    AddWaveformValueCalc<AddIsOnCalc<AddExposureDuration<AddStartTime<
        AddPulseDuration<AddPulseRepititionFrequency<Base<double> > > > > > >;

}

#endif  // include protector
