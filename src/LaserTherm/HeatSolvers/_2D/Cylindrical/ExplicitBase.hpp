#include <libField/Field.hpp>
#include <math.h>
#include <exception>
#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"
#include "LaserTherm/HeatSolvers/_2D/Cylindrical/FiniteDifferenceHeatSolver.hpp"

namespace BC = HeatSolvers::BoundaryConditions;
namespace FDHS = HeatSolvers::_2D::Cylindrical;

template <class IMP, class REAL>
class ExplicitBase : public FDHS::FiniteDifferenceHeatSolver<REAL> {
  public:
    // -------------------- PUBLIC CONSTRUCTORS --------------------
    // Constructors must be defined in derived classes
    // ---------------------- PUBLIC METHODS -----------------------
    void stepForward(REAL delta_t){
      REAL T1, T2, T3, T4, T5;
      REAL beta;
      for(int i = 0; i < zN; i++){
        for(int j = 0; j < rN; j++){
          // flags for boundary conditions
          bool rMax, rMin, zMax, zMin;
          rMax = j == rN-1;
          rMin = j == 0;
          zMax = i == zN-1;
          zMin = i == 0;
          // coefficient for calculations
          beta = delta_t / this->VHC[i][j];
          T1 = this->T[i][j] * static_cast<IMP*>(this)->A_n(i , j);
          if(rMin){
            // r=0 case
            // overwrite T1
            T1 = this->T[i][0]   * static_cast<IMP*>(this)->A_nR0(i , 0);
            T2 = this->T[i][1]   * static_cast<IMP*>(this)->B_nR0(i , 0);
            T3 = this->T[i][1]   * static_cast<IMP*>(this)->C_nR0(i , 0);
            //^^T[i][1] from symmetry about origin
          } else if(rMax){
            // r = Rmax
            switch(this->maxRBC.type){
              case BC::Type::Temperature:
                T2 = this->T[i][j-1] * static_cast<IMP*>(this)->B_n(i , rN-1);
                T3 = this->maxRBC.f  * static_cast<IMP*>(this)->C_n(i , rN-1);
                break;
              case BC::Type::HeatFlux:
                T2 = this->T[i][j-1] * static_cast<IMP*>(this)->B_n(i , rN-1);
                T3 = (this->T[i][rN-2] + static_cast<IMP*>(this)->get_dr(i, j) * this->maxRBC.f) * static_cast<IMP*>(this)->C_n(i , rN-1);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else {
            // if not on r boundary give 'normal T2 T3'
            T2 = this->T[i][j-1] * static_cast<IMP*>(this)->B_n(i , j);
            T3 = this->T[i][j+1] * static_cast<IMP*>(this)->C_n(i , j);
          }
          if(zMin){
            // z = 0
            switch(this->minZBC.type){
              case BC::Type::Temperature:
                // r = {0, rmax} caught by earlier ifs
                T4 = this->minZBC.f  * static_cast<IMP*>(this)->D_n(0 , j);
                T5 = this->T[1][j]   * static_cast<IMP*>(this)->E_n(0 , j);
                break;
              case BC::Type::HeatFlux:
                T4 = (this-> T[1][j] - static_cast<IMP*>(this)->get_dz(i,j) * this->minZBC.f)  * static_cast<IMP*>(this)->D_n(0 , j);
                T5 = this->T[1][j]   * static_cast<IMP*>(this)->E_n(0 , j);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else if(zMax){
            // z = Zmax
            switch(this->maxZBC.type){
              case BC::Type::Temperature:
                // r = {0, rmax} caught by earlier ifs
                T4 = this->T[i-1][j] * static_cast<IMP*>(this)->D_n(zN-1, j);
                T5 = this->maxZBC.f  * static_cast<IMP*>(this)->E_n(zN-1, j);
                break;
              case BC::Type::HeatFlux:
                T4 = this->T[i-1][j] * static_cast<IMP*>(this)->D_n(zN-1, j);
                T5 = (this->T[zN-2][j] + static_cast<IMP*>(this)->get_dz(i,j) * this->maxZBC.f)  * static_cast<IMP*>(this)->E_n(zN-1, j);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else {
            // if not on z boundary give 'normal T4 T5'
            T4 = this->T[i-1][j] * static_cast<IMP*>(this)->D_n(i , j);
            T5 = this->T[i+1][j] * static_cast<IMP*>(this)->E_n(i , j);
          }
          T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5 + this->A[i][j]);
          if(isnan(T_prime[i][j])){
            // make this a better error
            throw 42;
          }
        }
      }
      this->T += T_prime;
    }

    int get_rdims(){
      return this->rN;
    }

    int get_zdims(){
      return this->zN;
    }
  protected:
    Field<REAL, 2> T_prime;
    int rN;
    int zN;
  private:
    double percErr(REAL num, REAL correct){
      return abs((num - correct) / correct);
    }

    double randSampleErr(Field<REAL, 2>& f, REAL (*sol)(REAL, REAL, REAL), REAL (*err)(REAL, REAL), REAL t, int n){
      double e = 0.0;
      int i, j;
      REAL r, z;
      for(int i = 0; i < n; i++){
        i = rand() % f.size(0);
        j = rand() % f.size(1);
        r = f.getCoord(i, j)[1];
        z = f.getCoord(i, j)[1];
        e += err(f[i][j], sol(r, z, t));
      }
      return e;
    }

    double avgErr(Field<REAL, 2>& f, REAL (*sol)(REAL, REAL, REAL), REAL (*err)(REAL, REAL), REAL t){
      double e = 0.0;
      REAL r, z;
      for(int i = 0; i < f.size(0); i++){
        for(int j = 0; j < f.size(1); j++){
          r = f.getCoord(i, j)[1];
          z = f.getCoord(i, j)[0];
          e += err(f[i][j], sol(r, z, t));
        }
      }
      e /= f.size(0) * f.size(1) * 1.;
      return e;
    }

    double totalErr(Field<REAL, 2>& f, REAL (*sol)(REAL, REAL, REAL), REAL (*err)(REAL, REAL), REAL t){
      double e = 0.0;
      REAL r, z;
      for(int i = 0; i < f.size(0); i++){
        for(int j = 0; j < f.size(1); j++){
          r = f.getCoord(i, j)[1];
          z = f.getCoord(i, j)[0];
          e += err(f[i][j], sol(r, z, t));
        }
      }
      return e;
    }
};
