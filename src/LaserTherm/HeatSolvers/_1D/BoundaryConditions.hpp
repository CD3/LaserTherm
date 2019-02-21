#ifndef LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp
#define LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp

/** @file BoundaryConditions.hpp
  * @brief Boundary condition implementations.
  * @author C.D. Clark III
  * @date 02/15/19
  */

namespace HeatSolvers::_1D::BoundaryConditions {

enum class Type {Temperature, HeatFlux, None};

template<typename REAL>
struct FiniteDifference
{
  REAL f = 0;
  REAL dfdT = 0;
  Type type = Type::None;
};


template<typename REAL>
struct Sink
{
  template<typename BCType>
  void setBoundaryCondition( BCType& BC, const REAL& T, const REAL& t)
  {
    BC.type = Type::Temperature;
    BC.f = 0;
  }
};

template<typename REAL>
struct Insulator
{
  template<typename BCType>
  void setBoundaryCondition( BCType& BC, const REAL& T, const REAL& t)
  {
    BC.type = Type::HeatFlux;
    BC.f = 0;
    BC.dfdT = 0;
  }
};

template<typename REAL>
struct Convective
{
  REAL he;
  REAL Tinf;

  template<typename BCType>
  void setBoundaryCondition( BCType& BC, const REAL& T, const REAL& t)
  {
    BC.type = Type::HeatFlux;
    BC.f = he*(T - Tinf);
    BC.dfdT = he;
  }
};

struct Evaporative
{
};

struct Radiative
{
};


}



#endif // include protector
