#ifndef LaserTherm_HeatSources__1D_FiniteDifferenceHeatSource_hpp
#define LaserTherm_HeatSources__1D_FiniteDifferenceHeatSource_hpp

/** @file FiniteDifferenceHeatSource.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/24/19
  */

#include<Field.hpp>
#include<boost/signals2/signal.hpp>

namespace HeatSources::_1D {

template<typename REAL>
class FiniteDifferenceHeatSource
{
  public:
    Field<REAL,1> A;    ///< Source Term

    FiniteDifferenceHeatSource() = default;
    FiniteDifferenceHeatSource(size_t N)
    : A(N)
    {
    }

    // don't allow copy-construct or assign
    // user needs to call `reset` if they need to rconfigure size
    FiniteDifferenceHeatSource<REAL>(const FiniteDifferenceHeatSource<REAL>& other) = delete;
    FiniteDifferenceHeatSource<REAL>& operator-(const FiniteDifferenceHeatSource<REAL>& other) = delete;

    void reset(size_t N)
    {
      A.reset(N);
    }

};



}


#endif // include protector
