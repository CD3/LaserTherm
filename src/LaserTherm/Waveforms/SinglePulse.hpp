#ifndef LaserTherm_Waveforms_SinglePulse_hpp
#define LaserTherm_Waveforms_SinglePulse_hpp

/** @file SinglePulse.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/28/19
 */

namespace Waveforms
{
template<typename REAL>
using SinglePulse = AddWaveformValueCalc<
    AddIsOnCalc<AddStartTime<AddPulseDuration<Base<double> > > > >;

}

#endif  // include protector
