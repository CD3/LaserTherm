#ifndef LaserTherm_HeatSolvers__1D_FiniteDifferenceHeatSolver_hpp
#define LaserTherm_HeatSolvers__1D_FiniteDifferenceHeatSolver_hpp

/** @file FiniteDifferenceHeatSolver.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/13/19
  */

#include<Field.hpp>
#include<boost/signals2/signal.hpp>

#include "../BoundaryConditions.hpp"

namespace HeatSolvers::_1D {

template<typename REAL>
class FiniteDifferenceHeatSolver
{
  public:
    using BC = BoundaryConditions::FiniteDifference<REAL>;

    boost::signals2::signal< void(Field<REAL,1>&,const REAL&)  > sig_askSourceTerm;
    boost::signals2::signal< void(BC&,const REAL&,const REAL&) > sig_askMinBoundaryCondition;
    boost::signals2::signal< void(BC&,const REAL&,const REAL&) > sig_askMaxBoundaryCondition;




  public:
    FiniteDifferenceHeatSolver() = default;
    FiniteDifferenceHeatSolver(size_t N)
    : T(N),
      A(T.getCoordinateSystemPtr()),
      VHC(T.getCoordinateSystemPtr()),
      k(T.getCoordinateSystemPtr())
    {
    }

    Field<REAL,1> T;   ///< Temperature
    Field<REAL,1> A;   ///< Source Term
    Field<REAL,1> VHC; ///< Volumetric Heat Capacity
    Field<REAL,1> k;   ///< Conductivity

    BC minBC;
    BC maxBC;

};

}


#endif // include protector
