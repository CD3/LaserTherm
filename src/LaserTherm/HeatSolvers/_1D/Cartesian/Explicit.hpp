#ifndef LaserTherm_HeatSolvers__1D_Explicit_hpp
#define LaserTherm_HeatSolvers__1D_Explicit_hpp

/** @file Explicit.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/17/19
 */

#include <Eigen/Dense>

#include "../../../Utils/FiniteDifference.hpp"
#include "../FiniteDifferenceSolver.hpp"

namespace HeatSolvers::_1D::Cartesian
{
template<typename REAL>
class Explcit : public FiniteDifferenceSolver<REAL>
{
 public:
  using FiniteDifferenceSolver<REAL>::FiniteDifferenceSolver;
  using FiniteDifferenceSolver<REAL>::T;
  using FiniteDifferenceSolver<REAL>::A;
  using FiniteDifferenceSolver<REAL>::VHC;
  using FiniteDifferenceSolver<REAL>::k;

  using FiniteDifferenceSolver<REAL>::minBC;
  using FiniteDifferenceSolver<REAL>::maxBC;

  void stepForward(const REAL& dt);
  REAL calcMaxTimeStep() const;
};

template<typename REAL>
REAL Explicit<REAL>::calcMaxTimeStep() const
{
}

template<typename REAL>
void Explicit<REAL>::stepForward(const REAL& dt)
{
}

}  // namespace HeatSolvers::_1D

#endif  // include protector
