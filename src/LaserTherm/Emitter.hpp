#ifndef LaserTherm_Emitter_hpp
#define LaserTherm_Emitter_hpp

/** @file Emitter.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/28/19
  */

/**
 * An "emitter" is a heat source combined with a waveform. It knows when it is "on", and
 * it knows how to compute a source term.
 */

#include <Field.hpp>

template<class HeatSourceType, class WaveformType>
class Emitter : public HeatSourceType, public WaveformType
{

  protected:
    template<typename T>
    struct get_real_type {};

    template<typename REAL>
    struct get_real_type<Field<REAL,1>> { using type = REAL; };

  public:
    Emitter() = default;
    Emitter(size_t N):HeatSourceType(N) {}
    
    bool source_term_computed = false;

    using source_field_type = decltype(HeatSourceType::A);
    using real_type = typename get_real_type<source_field_type>::type;

    void addSourceTerm( source_field_type& A, const real_type& t);

};

template<class HeatSourceType, class WaveformType>
void Emitter<HeatSourceType, WaveformType>::addSourceTerm( source_field_type& A, const real_type& t)
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



#endif // include protector
