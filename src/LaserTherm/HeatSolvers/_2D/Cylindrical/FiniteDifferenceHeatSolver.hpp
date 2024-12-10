#pragma once
#ifndef LaserTherm_HeatSolvers__2D_Cylindrical_FiniteDifferenceHeatSolver_hpp
#define LaserTherm_HeatSolvers__2D_Cylindrical_FiniteDifferenceHeatSolver_hpp

/** @file FiniteDifferenceHeatSolver.hpp
 * @brief
 * @author C.D. Clark III
 */

#include <boost/signals2/signal.hpp>

#include <libField/Field.hpp>

#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"

namespace HeatSolvers::_2D::Cylindrical {
template <typename REAL> class FiniteDifferenceHeatSolver {
public:
  using BC = BoundaryConditions::FiniteDifference<REAL>;

  boost::signals2::signal<void(Field<REAL, 1> &, const REAL &)>
      sig_askSourceTerm;
  boost::signals2::signal<void(BC &, const REAL &, const REAL &)>
      sig_askRMinBoundaryCondition;
  boost::signals2::signal<void(BC &, const REAL &, const REAL &)>
      sig_askRMaxBoundaryCondition;
  boost::signals2::signal<void(BC &, const REAL &, const REAL &)>
      sig_askZMinBoundaryCondition;
  boost::signals2::signal<void(BC &, const REAL &, const REAL &)>
      sig_askZMaxBoundaryCondition;

public:
  FiniteDifferenceHeatSolver() = default;
  FiniteDifferenceHeatSolver(size_t Nz, size_t Nr)
      : T(Nz, Nr), A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()), k(T.getCoordinateSystemPtr()) {}

  FiniteDifferenceHeatSolver(const FiniteDifferenceHeatSolver &that)
      : T(that.T.size()), A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()), k(T.getCoordinateSystemPtr()) {
    for (int i = 0; i < this->T.size(); ++i) {
      this->T(i) = that.T(i);
      this->A(i) = that.A(i);
      this->VHC(i) = that.VHC(i);
      this->k(i) = that.k(i);
    }
  }

  auto
  operator=(FiniteDifferenceHeatSolver that) -> FiniteDifferenceHeatSolver {
    using std::swap;
    swap(T, that.T);
    swap(A, that.A);
    swap(VHC, that.VHC);
    swap(k, that.k);

    return *this;
  }

  // IMPORTANT
  // The order here matters!! Member variables initialized
  // in order of declaration. The constructor is passing
  // T.getCoordinateSystemPtr() to the constructors of A, VHC, and k. This means
  // the T object needs to be initialized _before_ the others.
  Field<REAL, 2> T;   ///< Temperature
  Field<REAL, 2> A;   ///< Source Term
  Field<REAL, 2> VHC; ///< Volumetric Heat Capacity
  Field<REAL, 2> k;   ///< Conductivity

  BC minZBC;
  BC maxZBC;
  BC minRBC;
  BC maxRBC;
};

} // namespace HeatSolvers::_2D::Cylindrical

#endif // include protector
