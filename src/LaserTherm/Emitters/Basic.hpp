#ifndef LaserTherm_Emitters_Basic_hpp
#define LaserTherm_Emitters_Basic_hpp


/** @file Basic.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/02/19
  */

/**
 * An "emitter" is a heat source combined with a waveform. It knows when it is "on", and
 * it knows how to compute a source term.
 */

#include <libField/Field.hpp>

#include "../Utils/TypeTraits.hpp"

namespace Emitters {

template<class HeatSourceType, class WaveformType>
class Basic: public HeatSourceType, public WaveformType
{

  public:
    Basic() = default;
    Basic(size_t N):HeatSourceType(N) {}
    
    bool source_term_computed = false;

    using source_field_type = decltype(HeatSourceType::A);
    using real_type = typename get_real_type<source_field_type>::type;

    void addSourceTerm( source_field_type& A, const real_type& t);
};

template<class HeatSourceType, class WaveformType>
void Basic<HeatSourceType, WaveformType>::addSourceTerm( source_field_type& A, const real_type& t)
{
  if( !this->isOn(t) )
    return;

  if( !source_term_computed )
  {
    HeatSourceType::compute();
    source_term_computed = true;
  }

  #pragma omp parallel for
  for(int i = 0; i < A.size(0); ++i)
  {
    A(i) += this->getWaveformValue(t)*HeatSourceType::A(i);
  }

}

}







#endif // include protector
