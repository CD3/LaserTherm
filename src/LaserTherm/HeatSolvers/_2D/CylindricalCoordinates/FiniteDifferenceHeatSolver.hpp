#pragma once

/** @file FiniteDifferenceHeatSolver.hpp
 * @brief
 * @author C.D. Clark III
 */

#include <boost/signals2/signal.hpp>

#include <libField/Field.hpp>

#include "../BoundaryConditions.hpp"

namespace HeatSolvers::_2D::CylindricalCoordinates
{
template<typename REAL>
class FiniteDifferenceHeatSolver
{
 public:
  using BC = BoundaryConditions::FiniteDifference<REAL>;

  boost::signals2::signal<void(Field<REAL, 1>&, const REAL&)> sig_askSourceTerm;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askRMinBoundaryCondition;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askRMaxBoundaryCondition;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askZMinBoundaryCondition;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askZMaxBoundaryCondition;

 public:
  FiniteDifferenceHeatSolver() = default;
  FiniteDifferenceHeatSolver(size_t N)
      : T(N),
        A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()),
        k(T.getCoordinateSystemPtr())
  {
  }

  FiniteDifferenceHeatSolver(const FiniteDifferenceHeatSolver& that)
      : T(that.T.size()),
        A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()),
        k(T.getCoordinateSystemPtr())
  {
    for (int i = 0; i < this->T.size(); ++i) {
      this->T(i)   = that.T(i);
      this->A(i)   = that.A(i);
      this->VHC(i) = that.VHC(i);
      this->k(i)   = that.k(i);
    }
  }

  FiniteDifferenceHeatSolver operator=(FiniteDifferenceHeatSolver that)
  {
    using std::swap;
    swap(T, that.T);
    swap(A, that.A);
    swap(VHC, that.VHC);
    swap(k, that.k);

    return *this;
  }

  Field<REAL, 1> T;    ///< Temperature
  Field<REAL, 1> A;    ///< Source Term
  Field<REAL, 1> VHC;  ///< Volumetric Heat Capacity
  Field<REAL, 1> k;    ///< Conductivity

  BC minBC;
  BC maxBC;
};

}  // namespace HeatSolvers::_1D

#endif  // include protector
