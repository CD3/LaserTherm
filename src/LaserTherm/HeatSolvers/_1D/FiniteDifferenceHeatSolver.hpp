#ifndef LaserTherm_HeatSolvers__1D_FiniteDifferenceHeatSolver_hpp
#define LaserTherm_HeatSolvers__1D_FiniteDifferenceHeatSolver_hpp

/** @file FiniteDifferenceHeatSolver.hpp
 * @brief
 * @author C.D. Clark III
 * @date 02/13/19
 */

#include <boost/signals2/signal.hpp>

#include <Field.hpp>

#include "../BoundaryConditions.hpp"

namespace HeatSolvers::_1D
{
template<typename REAL>
class FiniteDifferenceHeatSolver
{
 public:
  using BC = BoundaryConditions::FiniteDifference<REAL>;

  boost::signals2::signal<void(Field<REAL, 1>&, const REAL&)> sig_askSourceTerm;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askMinBoundaryCondition;
  boost::signals2::signal<void(BC&, const REAL&, const REAL&)>
      sig_askMaxBoundaryCondition;

 public:
  FiniteDifferenceHeatSolver() = default;
  FiniteDifferenceHeatSolver(size_t N)
      : T(N),
        A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()),
        k(T.getCoordinateSystemPtr())
  { }

  FiniteDifferenceHeatSolver(const FiniteDifferenceHeatSolver& that)
  :T(that.T.size()),
        A(T.getCoordinateSystemPtr()),
        VHC(T.getCoordinateSystemPtr()),
        k(T.getCoordinateSystemPtr())
  {
    for(int i = 0; i < this->T.size(); ++i)
    {
      this->T(i) = that.T(i);
      this->A(i) = that.A(i);
      this->VHC(i) = that.VHC(i);
      this->k(i) = that.k(i);
    }
  }

  FiniteDifferenceHeatSolver operator=(FiniteDifferenceHeatSolver that)
  {
    using std::swap;
    swap(T,that.T);
    swap(A,that.A);
    swap(VHC,that.VHC);
    swap(k,that.k);

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
