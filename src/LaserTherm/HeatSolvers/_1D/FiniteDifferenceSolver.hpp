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

    boost::signals2::signal< void(Field<REAL,1>&) > sig_askInitialTemperature;
    boost::signals2::signal< void(Field<REAL,1>&) > sig_askVolumetricHeatCapacity;
    boost::signals2::signal< void(Field<REAL,1>&) > sig_askConductivity;
    boost::signals2::signal< void(Field<REAL,1>&) > sig_askSourceTerm;
    //boost::signals2::signal< void(&) > sig_askBoundaryConditions;
    
    struct BC
    {
      REAL f = 0;
      REAL dfdT = 0;
      BoundaryConditions::Type type = BoundaryConditions::Type::None;
    };



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
