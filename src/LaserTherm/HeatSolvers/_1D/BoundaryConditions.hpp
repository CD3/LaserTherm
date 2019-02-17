#ifndef LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp
#define LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp

/** @file BoundaryConditions.hpp
  * @brief Boundary condition implementations.
  * @author C.D. Clark III
  * @date 02/15/19
  */

namespace HeatSolvers::_1D::BoundaryConditions {

enum class Type {Dirichlet, Neumann, None};

template<typename REAL>
struct FiniteDifference
{
  REAL f = 0;
  REAL dfdT = 0;
  Type type = Type::None;
};


struct Sink
{
  template<typename BCType, typename REAL, typename SolverType>
  void setBoundaryConditions( BCType& BC, REAL t, const SolverType& solver )
  {
    BC.type = Type::Dirichlet;
    BC.f = 0;
  }
};

struct Insulator
{
};

struct Convective
{
};

struct Evaporative
{
};

struct Radiative
{
};


}



#endif // include protector
