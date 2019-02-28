#ifndef LaserTherm_Waveforms_Mixins_hpp
#define LaserTherm_Waveforms_Mixins_hpp

/** @file Mixins.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/27/19
  */

#include <boost/hana.hpp>
#include <cmath>

namespace Waveforms {


template<typename REAL>
class Base
{
  public:
    using real_type = REAL;
};

#define MAKE_ADD_MEMBER_MIXIN( NAME ) \
template<typename BASE> \
class Add##NAME : public BASE \
{ \
  public: \
    using real_type = typename BASE::real_type; \
  protected: \
    real_type NAME; \
  public: \
    void set##NAME(const real_type& a){ NAME = a; } \
    real_type get##NAME() const { return NAME; } \
};\
template<typename T> \
constexpr auto has##NAME(T x) -> decltype(x.get##NAME(), std::true_type{}) {return{};} \
constexpr auto has##NAME(...) -> std::false_type {return{};}

MAKE_ADD_MEMBER_MIXIN( ExposureDuration );
MAKE_ADD_MEMBER_MIXIN( StartTime );
MAKE_ADD_MEMBER_MIXIN( PulseDuration );
MAKE_ADD_MEMBER_MIXIN( PulseRepititionFrequency );
MAKE_ADD_MEMBER_MIXIN( PulseNumber );

#undef MAKE_ADD_MEMBER_MIXIN



template<typename BASE>
class AddWaveformValueCalc: public BASE
{
  public:
    using real_type = typename BASE::real_type;
    real_type getWaveformValue( const real_type& time ) const
    {
      if( this->isOn(time) )
        return static_cast<real_type>(1);

      return static_cast<real_type>(0);
    }
};

template<typename BASE>
class AddIsOnCalc: public BASE
{
  public:
    using real_type = typename BASE::real_type;
    bool isOn( const real_type& time ) const
    {
      bool on = false;

      // if a start time is defined and the time passed in
      // is less than it, return false
      real_type startTime = 0;
      if constexpr( hasStartTime(BASE{}) )
        startTime = this->getStartTime();
      else
        on = true;

      if(time < startTime)
        return false;
      else
        on = true;


      // if an end time is definedand the time passed in
      // is greater than it, return false
      real_type endTime = -1;
      if constexpr( hasExposureDuration(BASE{}) )
        endTime = startTime + this->getExposureDuration();

      if( endTime > 0 && time > endTime )
        return false;




      if constexpr( hasPulseDuration(BASE{})
                && !hasPulseRepititionFrequency(BASE{}))
        endTime = startTime + this->getPulseDuration();

      if( endTime > 0 && time > endTime )
        return false;


      if constexpr( hasPulseDuration(BASE{})
                 && hasPulseRepititionFrequency(BASE{}))
      {
        real_type trash;
        real_type relative_pulse_duration = this->getPulseDuration() * this->getPulseRepititionFrequency();
        real_type relative_time = modf( time * this->getPulseRepititionFrequency(), &trash );
        if( relative_time > relative_pulse_duration )
          return false;
        else
          on = true;
      }


      return on;
    }
};



}



#endif // include protector
