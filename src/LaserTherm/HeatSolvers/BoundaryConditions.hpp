#ifndef LaserTherm_HeatSolvers_BoundaryConditions_hpp
#define LaserTherm_HeatSolvers_BoundaryConditions_hpp

/** @file BoundaryConditions.hpp
 * @brief Boundary condition implementations.
 * @author C.D. Clark III
 * @date 02/15/19
 */

namespace HeatSolvers::BoundaryConditions
{
  enum class Type { Temperature, HeatFlux, None };

  template<typename REAL>
    struct FiniteDifference {
      REAL f    = 0;
      REAL dfdT = 0;
      Type type = Type::None;
    };

  template<typename REAL>
    struct ConstantTemperature {
      REAL BoundaryTemperature = 0;
      ConstantTemperature(REAL aT) : BoundaryTemperature(aT) {}
      ConstantTemperature() = default;

      template<typename BCType>
        void setBoundaryCondition(BCType& BC, const REAL& aT = 0, const REAL& at = 0)
        {
          BC.type = Type::Temperature;
          BC.f    = BoundaryTemperature;
        }
    };

  template<typename REAL>
    using Sink = ConstantTemperature<REAL>;

  template<typename REAL>
    struct Insulator {
      Insulator() = default;

      template<typename BCType>
        void setBoundaryCondition(BCType& BC, const REAL& T, const REAL& t)
        {
          BC.type = Type::HeatFlux;
          BC.f    = 0;
          BC.dfdT = 0;
        }
    };

  template<typename REAL>
    struct Convective {
      REAL Tinf;
      REAL he;

      Convective(REAL aTinf, REAL ahe) : Tinf(aTinf), he(ahe) {}

      template<typename BCType>
        void setBoundaryCondition(BCType& BC, const REAL& T, const REAL& t)
        {
          BC.type = Type::HeatFlux;
          BC.f    = -he * (T - Tinf);
          BC.dfdT = -he;
        }
    };

  template<typename REAL>
    struct ConstantHeatFlux {
      REAL BoundaryHeatFlux;
      ConstantHeatFlux(REAL Q) : BoundaryHeatFlux(Q) {}
      ConstantHeatFlux() = default;

      template<typename BCType>
        void setBoundaryCondition(BCType& BC, const REAL& aT = 0, const REAL& at = 0)
        {
          BC.type = Type::HeatFlux;
          BC.f    = BoundaryHeatFlux;
          BC.dfdT = 0;
        }
    };

  struct Radiative {
  };

}  // namespace HeatSolvers::BoundaryConditions

#endif  // include protector
