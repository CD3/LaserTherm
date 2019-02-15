
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
  //public:

    using VectorType = Eigen::Matrix<REAL,Eigen::Dynamic,1>;

    VectorType Asub;
    VectorType Adiag;
    VectorType Asup;
    VectorType x;
    VectorType b;

  public:


    using FiniteDifferenceSolver<REAL>::FiniteDifferenceSolver;
    using FiniteDifferenceSolver<REAL>::T;
    using FiniteDifferenceSolver<REAL>::A;
    using FiniteDifferenceSolver<REAL>::VHC;
    using FiniteDifferenceSolver<REAL>::k;

    CrankNicholson(size_t N)
    :FiniteDifferenceSolver<REAL>(N),
     Asub(N),
     Adiag(N),
     Asup(N),
     x(N),
     b(N)
    { }

    void stepForward( const REAL& dt );



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
    b(i) = aR(i)     * T(i-1) 
         + bR(i,dt) * T(i) 
         + cR(i)     * T(i+1) ;
  }

  // LAST ELEMENT
  b(N-1) = aR(N-1)     * T(N-2) 
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
//

  // ____  _   _ ____  
  //|  _ \| | | / ___| 
  //| |_) | |_| \___ \ 
  //|  _ <|  _  |___) |
  //|_| \_\_| |_|____/ 
                   
  // MIN
  //if(minBCs != nullptr && minBCs->size() > 0 )
  //{
        //for( typename bcType::iterator it = minBCs->begin(); it != minBCs->end(); it++)
        //{
          //auto BC = it->second;
          //if( BC->get_type() == Neumann)
          //{
            //b(0)   += bRp(0, BC) * sol(0)
                   //+  cRp(0, BC) * sol(1)
                   //+   dp(0, BC);
          //}
     //else if( BC->get_type() == Dirichlet)
          //{
            //b(0)   += dp(0, BC);
          //}
        //}
  //}

  // MAX
  //if(maxBCs != nullptr && maxBCs->size() > 0 )
  //{
      //for( typename bcType::iterator it = maxBCs->begin(); it != maxBCs->end(); it++)
      //{
        //auto BC = it->second;
        //if( BC->get_type() == Neumann)
        //{
            //b(N-1)   += aRp(N-1, BC) * sol(N-2)
                     //+  bRp(N-1, BC) * sol(N-1)
                     //+   dp(N-1, BC);
        //}
   //else if( BC->get_type() == Dirichlet)
        //{
            //b(N-1)   += dp(N-1, BC);
        //}
      //}
  //}


  // _     _   _ ____  
  //| |   | | | / ___| 
  //| |   | |_| \___ \ 
  //| |___|  _  |___) |
  //|_____|_| |_|____/ 
                   
  // MIN
  //if(minBCs != nullptr && minBCs->size() > 0 )
  //{
    //for( typename bcType::iterator it = minBCs->begin(); it != minBCs->end(); it++)
    //{
      //auto BC = it->second;
      //if( BC->get_type() != Dirichlet)
      //{
        //Adiag(0)   += bLp(0, BC);
        //Asup(0)    += cLp(0, BC);
      //}
    //}
  //}

  // MAX
  //if(maxBCs != nullptr && maxBCs->size() > 0 )
  //{
    //for( typename bcType::iterator it = maxBCs->begin(); it != maxBCs->end(); it++)
    //{
      //auto BC = it->second;
      //if( BC->get_type() != Dirichlet)
      //{
        //Asub(N-1)  += aLp(N-1, BC);
        //Adiag(N-1) += bLp(N-1, BC);
      //}
    //}
  //}




//                                _                      
// ___  ___  _   _ _ __ ___ ___  | |_ ___ _ __ _ __ ___  
/// __|/ _ \| | | | '__/ __/ _ \ | __/ _ \ '__| '_ ` _ \ 
//\__ \ (_) | |_| | | | (_|  __/ | ||  __/ |  | | | | | |
//|___/\___/ \__,_|_|  \___\___|  \__\___|_|  |_| |_| |_|
                                                       

  //ADD PERFUSION TERM AND SOURCE TERM
  //if(this->sourceTermManifold != nullptr )
    //this->sourceTermManifold->set(0.0);
  //this->sig_updateSourceTerm( this->time, this->sourceTermManifold, &state );

  //#pragma omp parallel for 
  //for( int i = 0; i < N; i++ )
  //{
    //if( mat(i)->get_bloodFlowRate() > 0.0 )
      //b(i) += -2.0 * sol(i) * mat(i)->get_bloodFlowRate() * mat(i)->get_specificHeat();

    //// factor of 2 comes form Crank-Nicholson method
    //if(_sourceTerm != nullptr )
      //b(i) += -2.0 * (*_sourceTerm)(i);

    //// factor of 2 comes form Crank-Nicholson method
    //if(this->sourceTermManifold != nullptr )
      //b(i) += -2.0 * (*this->sourceTermManifold)(i);
  //}



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
REAL   CrankNicholson<REAL>::aLp( int _i )
{
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bLp( int _i )
{
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cLp( int _i )
{
}

template<typename REAL>
REAL   CrankNicholson<REAL>::aR(  int _i )
{

  REAL dzm,dzc,dzp;   // backward, center, and forward differences
  REAL dkdz;          // derivative of kappa (thermal conductivity )
  REAL km,kc,kp;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  kc = k(_i);
  km =  _i == 0           ? kc : k(_i-1);
  kp =  _i == k.size(0)-1 ? kc : k(_i+1);
    
  dkdz = firstDeriv_centFDmCoeff( k.getAxis(0), _i) * km
       + firstDeriv_centFDcCoeff( k.getAxis(0), _i) * kc
       + firstDeriv_centFDpCoeff( k.getAxis(0), _i) * kp;

 return (2*k(_i)- dkdz*dzp) / (dzc*dzm);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::bR(  int _i, REAL dt)
{

  REAL dzm,dzc,dzp;
  REAL dkdz;
  REAL km,kc,kp;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  kc = k(_i);
  km =  _i == 0           ? kc : k(_i-1);
  kp =  _i == k.size(0)-1 ? kc : k(_i+1);
    
  dkdz = firstDeriv_centFDmCoeff( k.getAxis(0), _i) * km
       + firstDeriv_centFDcCoeff( k.getAxis(0), _i) * kc
       + firstDeriv_centFDpCoeff( k.getAxis(0), _i) * kp;

 return beta(_i, dt) + (dkdz*(dzp - dzm) - 2*k(_i)) / (dzm*dzp);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::cR(  int _i )
{
  REAL dzm,dzc,dzp;   // backward, center, and forward differences
  REAL dkdz;          // derivative of kappa (thermal conductivity )
  REAL km,kc,kp;

  dzm = backDiff( k.getAxis(0), _i );
  dzc = centDiff( k.getAxis(0), _i );
  dzp = forwDiff( k.getAxis(0), _i );

  kc = k(_i);
  km =  _i == 0           ? kc : k(_i-1);
  kp =  _i == k.size(0)-1 ? kc : k(_i+1);
    
  dkdz = firstDeriv_centFDmCoeff( k.getAxis(0), _i) * km
       + firstDeriv_centFDcCoeff( k.getAxis(0), _i) * kc
       + firstDeriv_centFDpCoeff( k.getAxis(0), _i) * kp;

 return (2*k(_i)+ dkdz*dzm) / (dzc*dzp);
}

template<typename REAL>
REAL   CrankNicholson<REAL>::aRp( int _i )
{}

template<typename REAL>
REAL   CrankNicholson<REAL>::bRp( int _i )
{}

template<typename REAL>
REAL   CrankNicholson<REAL>::cRp( int _i )
{}

template<typename REAL>
REAL    CrankNicholson<REAL>::dp( int _i )
{}



}

#endif // include protector
