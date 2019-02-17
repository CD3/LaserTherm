
#ifndef LaserTherm_HeatSolvers__1D_CrankNicholson_hpp
#define LaserTherm_HeatSolvers__1D_CrankNicholson_hpp

/** @file CrankNicholson.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/13/19
  */

#include <Eigen/Dense>
#include "./FiniteDifferenceSolver.hpp"
#include "../../Utils/FiniteDifference.hpp"
#include "../../Utils/TriDiagonalSolver.hpp"


namespace HeatSolvers::_1D {

template<typename REAL>
class CrankNicholson : public FiniteDifferenceSolver<REAL>
{
  protected:

    using VectorType = Eigen::Matrix<REAL,Eigen::Dynamic,1>;

    VectorType Asub;
    VectorType Adiag;
    VectorType Asup;
    VectorType x;
    VectorType b;

  public:

    using FiniteDifferenceSolver<REAL>::sig_askSourceTerm;
    using FiniteDifferenceSolver<REAL>::sig_askMinBoundaryCondition;
    using FiniteDifferenceSolver<REAL>::sig_askMaxBoundaryCondition;

    using FiniteDifferenceSolver<REAL>::askSourceTerm;
    using FiniteDifferenceSolver<REAL>::askMinBoundaryCondition;
    using FiniteDifferenceSolver<REAL>::askMaxBoundaryCondition;

    using FiniteDifferenceSolver<REAL>::FiniteDifferenceSolver;
    using FiniteDifferenceSolver<REAL>::T;
    using FiniteDifferenceSolver<REAL>::A;
    using FiniteDifferenceSolver<REAL>::VHC;
    using FiniteDifferenceSolver<REAL>::k;

    using FiniteDifferenceSolver<REAL>::minBC;
    using FiniteDifferenceSolver<REAL>::maxBC;

    CrankNicholson(size_t N)
    :FiniteDifferenceSolver<REAL>(N),
     Asub(N),
     Adiag(N),
     Asup(N),
     x(N),
     b(N)
    { }

    void stepForward( const REAL& dt );
    REAL calcMaxTimeStep() const;




    REAL beta(  int _i, REAL dt);
    REAL   aL(  int _i );
    REAL   bL(  int _i, REAL dt);
    REAL   cL(  int _i );
    REAL   aLp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL   bLp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL   cLp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL   aR(  int _i );
    REAL   bR(  int _i, REAL dt);
    REAL   cR(  int _i );
    REAL   aRp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL   bRp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL   cRp( int _i ); // BoundaryConditionInterface<REAL>* _BC );
    REAL    dp( int _i ); // BoundaryConditionInterface<REAL>* _BC );

};


template<typename REAL>
REAL CrankNicholson<REAL>::calcMaxTimeStep( ) const
{
  // C-N is unconditionally stable, but we can get oscillations if
  // dt * (k/VHC) / dx = dt * k / (VHC * dx) > 1/2
  //
  // so, we shouldn't take time steps larger than
  //
  // dt = dx * VHZ / (2 * k)
  REAL dt = forwDiff( k.getAxis(0), 0 ) * VHC(0) / 2 / k(0);
  for( int i = 1; i < VHC.size(0); ++i)
  {
    REAL dt2 = forwDiff( k.getAxis(0), i ) * VHC(i) / 2 / k(i);
    if( dt2 < dt )
      dt = dt2;
  }

  return dt;
}

template<typename REAL>
void CrankNicholson<REAL>::stepForward( const REAL& dt )
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
  b(0) = bR(0,dt) * T(0)
       + cR(0)    * T(1);

  // INTERIOR ELEMENTS
  #pragma omp parallel for 
  for( int i = 1; i < N-1; i++)
  {
    b(i) = aR(i)    * T(i-1) 
         + bR(i,dt) * T(i) 
         + cR(i)    * T(i+1) ;
  }

  // LAST ELEMENT
  b(N-1) = aR(N-1)    * T(N-2) 
         + bR(N-1,dt) * T(N-1);
//        _    
//       / \   
//      / _ \  
//     / ___ \ 
//    /_/   \_\

  #pragma omp parallel for 
  for( int i = 0; i < N; ++i )
  {
    Asub(i)  = aL(i);
    Adiag(i) = bL(i, dt);
    Asup(i)  = cL(i);
  }


// _                           _                                        _ _ _   _                 
//| |__   ___  _   _ _ __   __| | __ _ _ __ _   _    ___ ___  _ __   __| (_) |_(_) ___  _ __  ___ 
//| '_ \ / _ \| | | | '_ \ / _` |/ _` | '__| | | |  / __/ _ \| '_ \ / _` | | __| |/ _ \| '_ \/ __|
//| |_) | (_) | |_| | | | | (_| | (_| | |  | |_| | | (_| (_) | | | | (_| | | |_| | (_) | | | \__ \
//|_.__/ \___/ \__,_|_| |_|\__,_|\__,_|_|   \__, |  \___\___/|_| |_|\__,_|_|\__|_|\___/|_| |_|___/
//                                          |___/                                                 

  // ____  _   _ ____  
  //|  _ \| | | / ___| 
  //| |_) | |_| \___ \ 
  //|  _ <|  _  |___) |
  //|_| \_\_| |_|____/ 
                   
  // TODO: we should actually get the BCs at 0 and dt if they are
  // time depdentnet. all of the "left" coefficients (aL, bL, cL) should
  // be multiplied by the BC at dt.
  sig_askMinBoundaryCondition( minBC, T(0), 0 );
  sig_askMaxBoundaryCondition( maxBC, T(N-1), 0 );

  // MIN
  if( minBC.type == BoundaryConditions::Type::Neumann)
    b(0) += bRp(0)*T(0) + cRp(0)*T(1) + dp(0);
  if( minBC.type == BoundaryConditions::Type::Dirichlet)
    b(0) += dp(0);

  // MAX
  if( maxBC.type == BoundaryConditions::Type::Neumann)
    b(N-1) += aRp(N-1)*T(N-2) + bRp(N-1)*T(N-1) + dp(N-1);
  if( maxBC.type == BoundaryConditions::Type::Dirichlet)
    b(N-1) += dp(N-1);
  


  // _     _   _ ____  
  //| |   | | | / ___| 
  //| |   | |_| \___ \ 
  //| |___|  _  |___) |
  //|_____|_| |_|____/ 

                   
  // MIN
  if( minBC.type == BoundaryConditions::Type::Neumann)
  {
    Adiag(0)   += bLp(0);
    Asup(0)    += cLp(0);
  }

  // MAX
  if( maxBC.type == BoundaryConditions::Type::Neumann)
  {
    Asub(N-1)  += aLp(N-1);
    Adiag(N-1) += bLp(N-1);
  }




//                                _                      
// ___  ___  _   _ _ __ ___ ___  | |_ ___ _ __ _ __ ___  
/// __|/ _ \| | | | '__/ __/ _ \ | __/ _ \ '__| '_ ` _ \ 
//\__ \ (_) | |_| | | | (_|  __/ | ||  __/ |  | | | | | |
//|___/\___/ \__,_|_|  \___\___|  \__\___|_|  |_| |_| |_|
                                                       

  // crank-nicholson averages the right-hand side of the heat
  // equation. the source term we use should be the sum
  // of current and next timetep source terms

  sig_askSourceTerm( A, 0 );
  #pragma omp parallel for 
  for( int i = 0; i < N; ++i )
    b(i) += A(i);

  sig_askSourceTerm( A, dt );
  #pragma omp parallel for 
  for( int i = 0; i < N; ++i )
    b(i) += A(i);



//           _           
// ___  ___ | |_   _____ 
/// __|/ _ \| \ \ / / _ \
//\__ \ (_) | |\ V /  __/
//|___/\___/|_| \_/ \___|
                       

  TriDiagonalSolver<Thomas>::Solve( Asub, Adiag, Asup, b, x );

  #pragma omp parallel for 
  for( int i = 0; i < N; i++ )
  {
    T(i) = x(i);
  }

  return;

}


template<typename REAL>
REAL CrankNicholson<REAL>::beta(  int _i, REAL dt)
{
  return (2. * VHC(_i) ) / dt;
}

template<typename REAL>
REAL   CrankNicholson<REAL>::aL(  int _i )
{
  return -aR( _i );
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bL(  int _i, REAL dt)
{
  return 2*beta(_i, dt) - bR( _i, dt );
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cL(  int _i )
{
  return -cR( _i );
}

template<typename REAL>
REAL   CrankNicholson<REAL>::aR(  int _i )
{

  REAL dzm,dzc,dzp;
  REAL dkdz;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  dkdz = firstDeriv_centFD(k,_i);

  return (2*k(_i)- dkdz*dzp) / (dzc*dzm);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bR(  int _i, REAL dt)
{
  REAL dzm,dzc,dzp;
  REAL dkdz;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  dkdz = firstDeriv_centFD(k,_i);

  return beta(_i, dt) + (dkdz*(dzp - dzm) - 2*k(_i)) / (dzm*dzp);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cR(  int _i )
{
  REAL dzm,dzc,dzp;
  REAL dkdz;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  dkdz = firstDeriv_centFD(k,_i);

  return (2*k(_i)+ dkdz*dzm) / (dzc*dzp);
}



template<typename REAL>
REAL   CrankNicholson<REAL>::aLp( int _i )
{
  return cL(_i);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bLp( int _i )
{
  return bRp( _i );
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cLp( int _i )
{
  return aL(_i);
}


template<typename REAL>
REAL   CrankNicholson<REAL>::aRp( int _i )
{
  return cR(_i);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bRp( int _i )
{
  if( _i == 0  )
    return -aL(_i) * centDiff( T.getAxis(0), _i ) * minBC.dfdT;
  else if( _i == T.size(0) - 1)
    return  cL(_i) * centDiff( T.getAxis(0), _i ) * maxBC.dfdT;

  return REAL(0);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cRp( int _i )
{
  return aR(_i);
}

template<typename REAL>
REAL    CrankNicholson<REAL>::dp( int _i )
{
  if( _i == 0 )
  {
    if( minBC.type == BoundaryConditions::Type::Neumann)
      return  (aL(_i) - aR(_i)) * centDiff( T.getAxis(0), _i ) * minBC.f;
    if( minBC.type == BoundaryConditions::Type::Dirichlet )
      return  (aR(_i) - aL(_i)) * minBC.f;
  }
  else if( _i == T.size(0)-1 )
  {
    if( maxBC.type == BoundaryConditions::Type::Neumann)
      return  (cL(_i) - cR(_i)) * centDiff( T.getAxis(0), _i ) * maxBC.f;
    if( maxBC.type == BoundaryConditions::Type::Dirichlet )
      return  (cR(_i) - cL(_i)) * maxBC.f;

  }

  return REAL(0);
}



}

#endif // include protector
