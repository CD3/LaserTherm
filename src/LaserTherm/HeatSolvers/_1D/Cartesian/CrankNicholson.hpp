
#ifndef LaserTherm_HeatSolvers__1D_CrankNicholson_hpp
#define LaserTherm_HeatSolvers__1D_CrankNicholson_hpp

/** @file CrankNicholson.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/13/19
 */

#include <Eigen/Dense>

#include "../../../Utils/FiniteDifference.hpp"
#include "../../../Utils/TriDiagonalSolver.hpp"
#include "../../../Utils/TypeTraits.hpp"
#include "../FiniteDifferenceHeatSolver.hpp"

namespace HeatSolvers::_1D::Cartesian
{
template<typename REAL>
class CrankNicholson : public FiniteDifferenceHeatSolver<REAL>
{
 protected:
  using VectorType = Eigen::Matrix<REAL, Eigen::Dynamic, 1>;

  VectorType Asub;
  VectorType Adiag;
  VectorType Asup;
  VectorType x;
  VectorType b;

 public:
  using FiniteDifferenceHeatSolver<REAL>::sig_askSourceTerm;
  using FiniteDifferenceHeatSolver<REAL>::sig_askMinBoundaryCondition;
  using FiniteDifferenceHeatSolver<REAL>::sig_askMaxBoundaryCondition;

  using FiniteDifferenceHeatSolver<REAL>::FiniteDifferenceHeatSolver;
  using FiniteDifferenceHeatSolver<REAL>::T;
  using FiniteDifferenceHeatSolver<REAL>::A;
  using FiniteDifferenceHeatSolver<REAL>::VHC;
  using FiniteDifferenceHeatSolver<REAL>::k;

  using BC = BoundaryConditions::FiniteDifference<REAL>;
  using FiniteDifferenceHeatSolver<REAL>::minBC;
  using FiniteDifferenceHeatSolver<REAL>::maxBC;

  CrankNicholson() = default;
  CrankNicholson(size_t N)
      : FiniteDifferenceHeatSolver<REAL>(N),
        Asub(N),
        Adiag(N),
        Asup(N),
        x(N),
        b(N)
  {
  }
  CrankNicholson(const CrankNicholson& that) = default;

  void stepForward(const REAL& dt);
  REAL calcMaxTimeStep() const;

  REAL beta(int _i, REAL dt);
  REAL aL(int _i);
  REAL bL(int _i, REAL dt);
  REAL cL(int _i);
  REAL aR(int _i);
  REAL bR(int _i, REAL dt);
  REAL cR(int _i);
  REAL aLp(int _i);
  REAL bLp(int _i);
  REAL cLp(int _i);
  REAL aRp(int _i);
  REAL bRp(int _i);
  REAL cRp(int _i);
  REAL dp1(int _i);
  REAL dp2(int _i);
};

template<typename REAL>
REAL CrankNicholson<REAL>::calcMaxTimeStep() const
{
  // C-N is unconditionally stable, but we can get oscillations if
  // dt * (k/VHC) / dx = dt * k / (VHC * dx) > 1/2
  //
  // so, we shouldn't take time steps larger than
  //
  // dt = dx * VHZ / (2 * k)
  REAL dt = forwDiff(k.getAxis(0), 0) * VHC(0) / 2 / k(0);
  for (int i = 1; i < VHC.size(0); ++i) {
    REAL dt2 = forwDiff(k.getAxis(0), i) * VHC(i) / 2 / k(i);
    if (dt2 < dt) dt = dt2;
  }

  return dt;
}

template<typename REAL>
void CrankNicholson<REAL>::stepForward(const REAL& dt)
{
  int N = T.size(0);

  ////////////////////////////////////
  ////////////////////////////////////
  ////////////////////////////////////
  //    _                    _      //
  //   / \   __  __  _____  | |__   //
  //  / _ \  \ \/ / |_____| | '_ \  //
  // / ___ \  >  <  |_____| | |_) | //
  ///_/   \_\/_/\_\         |_.__/  //
  ////////////////////////////////////
  ////////////////////////////////////
  ////////////////////////////////////

  //            _
  //           | |__
  //           | '_ \
//           | |_) |
  //           |_.__/

  // FIRST ELEMENT
  b(0) = bR(0, dt) * T(0) + cR(0) * T(1);

// INTERIOR ELEMENTS
#pragma omp parallel for
  for (int i = 1; i < N - 1; i++) {
    b(i) = aR(i) * T(i - 1) + bR(i, dt) * T(i) + cR(i) * T(i + 1);
  }

  // LAST ELEMENT
  b(N - 1) = aR(N - 1) * T(N - 2) + bR(N - 1, dt) * T(N - 1);
  //        _
  //       / \
//      / _ \
//     / ___ \
//    /_/   \_\

#pragma omp parallel for
  for (int i = 0; i < N; ++i) {
    Asub(i)  = aL(i);
    Adiag(i) = bL(i, dt);
    Asup(i)  = cL(i);
  }

  // _                           _                                        _ _ _
  // _
  //| |__   ___  _   _ _ __   __| | __ _ _ __ _   _    ___ ___  _ __   __| (_)
  //|_(_) ___  _ __  ___ | '_ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | |  / __/
  //_ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
  //| |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | | (_| (_) | | | | (_| | |
  //|_| | (_) | | | \__ \
//|_.__/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |  \___\___/|_|
  //|_|\__,_|_|\__|_|\___/|_| |_|___/
  //                                          |___/

  // ____  _   _ ____
  //|  _ \| | | / ___|
  //| |_) | |_| \___ \
  //|  _ <|  _  |___) |
  //|_| \_\_| |_|____/

  // the dp coefficients need to be evaluated at the current
  // time
  sig_askMinBoundaryCondition(minBC, T(0), 0);
  sig_askMaxBoundaryCondition(maxBC, T(N - 1), 0);

  b(0) += dp1(0);          // MIN
  b(N - 1) += dp1(N - 1);  // MAX

  // now we can evaluate the BC's at the next time
  sig_askMinBoundaryCondition(minBC, T(0), dt);
  sig_askMaxBoundaryCondition(maxBC, T(N - 1), dt);

  // MIN
  b(0) += dp2(0);
  if (minBC.type == BoundaryConditions::Type::HeatFlux)
    b(0) += bRp(0) * T(0) + cRp(0) * T(1);  //  + dp2(0)

  // MAX
  b(N - 1) += dp2(N - 1);
  if (maxBC.type == BoundaryConditions::Type::HeatFlux)
    b(N - 1) += aRp(N - 1) * T(N - 2) + bRp(N - 1) * T(N - 1);  // + dp2(N-1);

  // _     _   _ ____
  //| |   | | | / ___|
  //| |   | |_| \___ \
  //| |___|  _  |___) |
  //|_____|_| |_|____/

  // MIN
  if (minBC.type == BoundaryConditions::Type::HeatFlux) {
    Adiag(0) += bLp(0);
    Asup(0) += cLp(0);
  }

  // MAX
  if (maxBC.type == BoundaryConditions::Type::HeatFlux) {
    Asub(N - 1) += aLp(N - 1);
    Adiag(N - 1) += bLp(N - 1);
  }

  //                                _
  // ___  ___  _   _ _ __ ___ ___  | |_ ___ _ __ _ __ ___
  /// __|/ _ \| | | | '__/ __/ _ \ | __/ _ \ '__| '_ ` _ \
//\__ \ (_) | |_| | | | (_|  __/ | ||  __/ |  | | | | | |
  //|___/\___/ \__,_|_|  \___\___|  \__\___|_|  |_| |_| |_|

  // crank-nicholson averages the right-hand side of the heat
  // equation. the source term we use should be the sum
  // of current and next timetep source terms

  sig_askSourceTerm(A, 0);
#pragma omp parallel for
  for (int i = 0; i < N; ++i) b(i) += A(i);

  sig_askSourceTerm(A, dt);
#pragma omp parallel for
  for (int i = 0; i < N; ++i) b(i) += A(i);

  //           _
  // ___  ___ | |_   _____
  /// __|/ _ \| \ \ / / _ \
//\__ \ (_) | |\ V /  __/
  //|___/\___/|_| \_/ \___|

  TriDiagonalSolver<Thomas>::Solve(Asub, Adiag, Asup, b, x);

#pragma omp parallel for
  for (int i = 0; i < N; i++) T(i) = x(i);

  return;
}

template<typename REAL>
REAL CrankNicholson<REAL>::beta(int _i, REAL dt)
{
  return (2. * VHC(_i)) / dt;
}

template<typename REAL>
REAL CrankNicholson<REAL>::aL(int _i)
{
  return -aR(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::bL(int _i, REAL dt)
{
  return 2 * beta(_i, dt) - bR(_i, dt);
}

template<typename REAL>
REAL CrankNicholson<REAL>::cL(int _i)
{
  return -cR(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::aR(int _i)
{
  REAL dzm, dzc, dzp;
  REAL dkdz;

  dzm = backDiff(k.getAxis(0), _i);
  dzc = centDiff(k.getAxis(0), _i);
  dzp = forwDiff(k.getAxis(0), _i);

  dkdz = firstDeriv_centFD(k, _i);

  return (2 * k(_i) - dkdz * dzp) / (dzc * dzm);
}

template<typename REAL>
REAL CrankNicholson<REAL>::bR(int _i, REAL dt)
{
  REAL dzm, dzc, dzp;
  REAL dkdz;

  dzm = backDiff(k.getAxis(0), _i);
  dzc = centDiff(k.getAxis(0), _i);
  dzp = forwDiff(k.getAxis(0), _i);

  dkdz = firstDeriv_centFD(k, _i);

  return beta(_i, dt) + (dkdz * (dzp - dzm) - 2 * k(_i)) / (dzm * dzp);
}

template<typename REAL>
REAL CrankNicholson<REAL>::cR(int _i)
{
  REAL dzm, dzc, dzp;
  REAL dkdz;

  dzm = backDiff(k.getAxis(0), _i);
  dzc = centDiff(k.getAxis(0), _i);
  dzp = forwDiff(k.getAxis(0), _i);

  dkdz = firstDeriv_centFD(k, _i);

  return (2 * k(_i) + dkdz * dzm) / (dzc * dzp);
}

template<typename REAL>
REAL CrankNicholson<REAL>::aLp(int _i)
{
  return cL(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::bLp(int _i)
{
  return bRp(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::cLp(int _i)
{
  return aL(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::aRp(int _i)
{
  return cR(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::bRp(int _i)
{
  if (_i == 0)
    return -aL(_i) * centDiff(T.getAxis(0), _i) *
           (-minBC.dfdT / k(_i));  // flip sign for min boundary
  else if (_i == T.size(0) - 1)
    return cL(_i) * centDiff(T.getAxis(0), _i) * (maxBC.dfdT / k(_i));

  return REAL(0);
}

template<typename REAL>
REAL CrankNicholson<REAL>::cRp(int _i)
{
  return aR(_i);
}

template<typename REAL>
REAL CrankNicholson<REAL>::dp1(int _i)
{
  if (_i == 0) {
    if (minBC.type == BoundaryConditions::Type::HeatFlux)
      return -aR(_i) * centDiff(T.getAxis(0), _i) *
             (-minBC.f / k(_i));  // flip sign for min boundary
    if (minBC.type == BoundaryConditions::Type::Temperature)
      return aR(_i) * minBC.f;
  } else if (_i == T.size(0) - 1) {
    if (maxBC.type == BoundaryConditions::Type::HeatFlux)
      return cR(_i) * centDiff(T.getAxis(0), _i) * maxBC.f / k(_i);
    if (maxBC.type == BoundaryConditions::Type::Temperature)
      return cR(_i) * maxBC.f;
  }

  return REAL(0);
}

template<typename REAL>
REAL CrankNicholson<REAL>::dp2(int _i)
{
  if (_i == 0) {
    if (minBC.type == BoundaryConditions::Type::HeatFlux)
      return aL(_i) * centDiff(T.getAxis(0), _i) *
             (-minBC.f / k(_i));  // flip sign for min boundary
    if (minBC.type == BoundaryConditions::Type::Temperature)
      return -aL(_i) * minBC.f;
  } else if (_i == T.size(0) - 1) {
    if (maxBC.type == BoundaryConditions::Type::HeatFlux)
      return -cL(_i) * centDiff(T.getAxis(0), _i) * maxBC.f / k(_i);
    if (maxBC.type == BoundaryConditions::Type::Temperature)
      return -cL(_i) * maxBC.f;
  }

  return REAL(0);
}

}  // namespace HeatSolvers::_1D::Cartesian

template<typename REAL>
struct get_real_type<HeatSolvers::_1D::Cartesian::CrankNicholson<REAL> > {
  using type = REAL;
};

#endif  // include protector
