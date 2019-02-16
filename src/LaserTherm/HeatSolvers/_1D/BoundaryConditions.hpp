#ifndef LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp
#define LaserTherm_HeatSolvers__1D_BoundaryConditions_hpp

/** @file BoundaryConditions.hpp
  * @brief Boundary condition implementations.
  * @author C.D. Clark III
  * @date 02/15/19
  */

namespace HeatSolvers::_1D::BoundaryConditions {

enum class Type {Dirichlet, Neumann, None};

struct Sink
{
};

struct Insulator
{
};



}



#endif // include protector
