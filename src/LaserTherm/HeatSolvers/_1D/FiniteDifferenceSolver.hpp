#ifndef LaserTherm_HeatSolvers__1D_FiniteDifferenceSolver_hpp
#define LaserTherm_HeatSolvers__1D_FiniteDifferenceSolver_hpp

/** @file FiniteDifferenceSolver.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/13/19
  */

#include<Field.hpp>
#include<boost/signals2/signal.hpp>

#include "./BoundaryConditions.hpp"

namespace HeatSolvers::_1D {

template<typename REAL>
class FiniteDifferenceSolver
{
  public:
    using BC = BoundaryConditions::FiniteDifference<REAL>;

    boost::signals2::signal< void(Field<REAL,1>&,const REAL&)  > sig_askSourceTerm;
    boost::signals2::signal< void(BC&,const REAL&,const REAL&) > sig_askMinBoundaryCondition;
    boost::signals2::signal< void(BC&,const REAL&,const REAL&) > sig_askMaxBoundaryCondition;


    // provide functions to connect signals in case we want to add multiple signatures
    // or change the signal/slot backend at some point.
    template<typename F>
    auto askSourceTerm( F&& f ){ return sig_askSourceTerm.connect(std::forward(f)); }
    template<typename F>
    auto askMinBoundaryCondition( F&& f ){ return sig_askMinBoundaryCondition.connect(std::forward(f)); }
    template<typename F>
    auto askMaxBoundaryCondition( F&& f ){ return sig_askMaxBoundaryCondition.connect(std::forward(f)); }
    



  public:
    FiniteDifferenceSolver() = default;
    FiniteDifferenceSolver(size_t N)
    : T(N),
      A(T.getCoordinateSystemPtr()),
      VHC(T.getCoordinateSystemPtr()),
      k(T.getCoordinateSystemPtr())
    {
    }

    void initialize()
    {
      sig_askInitialTemperature( T );
      sig_askVolumetricHeatCapacity( VHC );
      sig_askConductivity( k );
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
