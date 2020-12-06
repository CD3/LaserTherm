#ifndef LaserTherm_HeatSources__1D_BeersLaw_hpp
#define LaserTherm_HeatSources__1D_BeersLaw_hpp

/** @file BeersLaw.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/24/19
 */

#include "./FiniteDifferenceHeatSource.hpp"

namespace HeatSources::_1D
{
template<typename REAL>
class BeersLaw : public FiniteDifferenceHeatSource<REAL>
{
 public:
  REAL           E0 = 1;  ///< Incident irradiance
  Field<REAL, 1> mu_a;    ///< Absorption coefficient
  using FiniteDifferenceHeatSource<REAL>::A;

  BeersLaw() = default;
  BeersLaw(size_t N)
      : FiniteDifferenceHeatSource<REAL>(N), mu_a(A.getCoordinateSystemPtr())
  {
  }

  void setIrradiance(REAL a_E0) { E0 = a_E0; };
  REAL getIrradiance() const { return E0; }

  void reset(size_t N)
  {
    FiniteDifferenceHeatSource<REAL>::reset(N);
    mu_a.reset(A.getCoordinateSystemPtr());
  }

  void compute();
};

template<typename REAL>
void BeersLaw<REAL>::compute()
{
  A(0)   = mu_a(0) * E0;
  REAL l = 0;
  for (int i = 1; i < A.size(); ++i) {
    l += mu_a(i) * (A.getAxis(0)[i] - A.getAxis(0)[i - 1]);
    A(i) = mu_a(i) * E0 * exp(-l);
  }
}

}  // namespace HeatSources::_1D

#endif  // include protector
