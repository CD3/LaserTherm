#ifndef LaserTherm_Waveforms_ContinuousWave_hpp
#define LaserTherm_Waveforms_ContinuousWave_hpp

/** @file ContinuousWave.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/28/19
  */

#include "./Mixins.hpp"

namespace Waveforms {

template<typename REAL>
using ContinuousWave = 
      AddWaveformValueCalc<
      AddIsOnCalc<
      AddExposureDuration<
      AddStartTime<
      Base<REAL> > > > > ;

}


#endif // include protector
