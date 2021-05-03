#include <libField/Field.hpp>
#include <math.h>
#include <exception>
#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"
#include "LaserTherm/HeatSolvers/_2D/Cylindrical/FiniteDifferenceHeatSolver.hpp"

namespace BC = HeatSolvers::BoundaryConditions;
namespace FDHS = HeatSolvers::_2D::Cylindrical;

template <class REAL>
class Explicit : public FDHS::FiniteDifferenceHeatSolver<REAL> {
  public:
    // -------------------- PUBLIC CONSTRUCTORS --------------------
    Explicit(size_t _zN, size_t _rN){
      this->rN = _rN;
      this->zN = _zN;
      this->T.reset(zN, rN);
      this->A.reset(this->T.getCoordinateSystemPtr());
      this->VHC.reset(this->T.getCoordinateSystemPtr());
      this->k.reset(this->T.getCoordinateSystemPtr());
      this->T_prime.reset(this->T.getCoordinateSystemPtr());
    }
    Explicit(){
      // Default constructor to be safe
    }
    // ---------------------- PUBLIC METHODS -----------------------
    // delta_t is a step forward in time in seconds
    void stepForward(REAL delta_t){
      REAL T1, T2, T3, T4, T5;
      REAL beta;
      for(int i = 0; i < zN; i++){
        for(int j = 0; j < rN; j++){
          beta = delta_t / this->VHC[i][j];
          T1 = this->T[i][j] * A_n(i , j);
          if(j == 0){
            // r=0 case
            // overwrite T1
            T1 = this->T[i][0]   * A_nR0(i , 0);
            T2 = this->T[i][1]   * B_nR0(i , 0);
            T3 = this->T[i][1]   * C_nR0(i , 0);
            //^^T[i][1] from symmetry about origin
          } else if(j == rN-1){
            // r = Rmax
            switch(this->maxRBC.type){
              case BC::Type::Temperature:
                T2 = this->maxRBC.f  * B_n(i , rN-1);
                T3 = this->T[i][j-1] * C_n(i , rN-1);
                break;
              case BC::Type::HeatFlux:
                T2 = (this->T[i][rN-2] + get_dr(i, j) * this->maxRBC.f)  * B_n(i , rN-1);
                T3 = this->T[i][j-1] * C_n(i , rN-1);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else {
            // if not on r boundary give 'normal T2 T3'
            T2 = this->T[i][j+1] * B_n(i , j);
            T3 = this->T[i][j-1] * C_n(i , j);
          }
          if(i == 0){
            // z = 0
            switch(this->minZBC.type){
              case BC::Type::Temperature:
                // r = {0, rmax} caught by earlier ifs
                T4 = this->minZBC.f  * D_n(0 , j);
                T5 = this->T[1][j]   * E_n(0 , j);
                break;
              case BC::Type::HeatFlux:
                T4 = (this-> T[1][j] - get_dz(i,j) * this->minZBC.f)  * D_n(0 , j);
                T5 = this->T[1][j]   * E_n(0 , j);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else if(i == zN-1){
            // z = Zmax
            switch(this->maxZBC.type){
              case BC::Type::Temperature:
                // r = {0, rmax} caught by earlier ifs
                T4 = this->T[i-1][j] * D_n(zN-1, j);
                T5 = this->maxZBC.f  * E_n(zN-1, j);
                break;
              case BC::Type::HeatFlux:
                T4 = this->T[i-1][j] * D_n(zN-1, j);
                T5 = (this->T[zN-2][j] + get_dz(i,j) * this->maxZBC.f)  * E_n(zN-1, j);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              default:
                throw 42;
            }
          } else {
            // if not on z boundary give 'normal T4 T5'
            T4 = this->T[i-1][j] * D_n(i , j);
            T5 = this->T[i+1][j] * E_n(i , j);
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
    // ----------------- PUBLIC MEMBER VARIABLES -------------------

  protected:
    // ---------------- PROTECTED MEMBER VARIABLES  ----------------
    // --------------------- PROTECTED METHODS ---------------------
    // get forward r spacing (backwards at rmax)
    REAL get_dr(int i, int j){
      REAL dr;
      if(j == rN - 1){
        dr = this->k.getCoord(i, j)[1] - this->k.getCoord(i, j - 1)[1];
      } else {
        dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      }
      return dr;
    }

    // get forward z spacing (backwards at zmax)
    REAL get_dz(int i, int j){
      REAL dz;
      if(i == zN - 1){
        dz = this->k.getCoord(i, j)[0] - this->k.getCoord(i - 1, j)[0];
      } else {
        dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      }
      return dz;
    }

    REAL forward_finite_r(Field<REAL, 2>& f, int i, int j){
      return f[i][j+1] - f[i][j] / get_dr(i, j);
    }
    REAL backward_finite_r(Field<REAL, 2>& f, int i, int j){
      return f[i][j] - f[i][j-1] / get_dr(i, j);
    }
    REAL central_finite_r(Field<REAL, 2>& f, int i, int j){
      return f[i][j+1] - f[i][j-1] / (2 * get_dr(i,j));
    }
    REAL forward_finite_z(Field<REAL, 2>& f, int i, int j){
      return f[i+1][j] - f[i][j] / get_dz(i, j);
    }
    REAL backward_finite_z(Field<REAL, 2>& f, int i, int j){
      return f[i][j] - f[i-1][j] / get_dz(i, j);
    }
    REAL central_finite_z(Field<REAL, 2>& f, int i, int j){
      return f[i+1][j] - f[i-1][j] / (2 * get_dz(i, j));
    }
    // Finite Difference deriv for kappa wrt r
    REAL dk_dr(int i, int j){
      if (j == 0){
        return forward_finite_r(this->k, i, j);
      }
      else if (j == rN-1){
        return backward_finite_r(this->k, i, j);
      }
      else{
        return central_finite_r(this->k, i, j);
      }
    }
    // Finite Difference deriv for kappa wrt z
    REAL dk_dz(int i, int j){
      if (i == 0){
        return forward_finite_z(this->k, i, j);
      }
      else if (i == zN-1){
        return backward_finite_z(this->k, i, j);
      }
      else{
        return central_finite_z(this->k, i, j);
      }
    }
    // Calculate Coefficent for T^n_(z, r)
    REAL A_n(int i, int j){
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      // -2 * kappa / dr**2
      REAL T1 = -2 * this->k[i][j] / pow(dr, 2);
      // -2 * kappa / dz**2
      REAL T2 = -2 * this->k[i][j] / pow(dz, 2);
      return T1 + T2;
    }
    // Calculate Coefficent for T^n_(z, r+1)
    REAL B_n(int i, int j){
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      // kappa / 2 r dr
      REAL T1 = this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      REAL T2 = this->k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      REAL T3 = dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }

    // Calculate Coefficent for T^n_(z, r-1)
    REAL C_n(int i, int j){
      //return (N) 0;
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      // - kappa / 2r * dr
      REAL T1 = - this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      REAL T2 = this->k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      REAL T3 = - dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }

    // Calculate Coefficent for T^n_(z-1, r)
    REAL D_n(int i, int j){
      //return (N) 0;
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      // kappa / dz**2
      REAL T1 = this->k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      REAL T2 = - dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }

    // Calculate Coefficent for T^n_(z+1, r)
    REAL E_n(int i, int j){
      //return (N) 0;
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      // kappa / dz**2
      REAL T1 = this->k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      REAL T2 = dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }

    // Calculate Coefficent for T^n_(z+1, 0)
    REAL A_nR0(int i, int j){
      REAL dr = get_dr(i, j);
      REAL dz = get_dz(i, j);
      REAL T1 = (-4 * this->k[i][j]) / (dr * dr);
      REAL T2 = (-2 * this->k[i][j]) / (dz * dz);
      return T1 + T2;
    }

    // Calculate Coefficent for T^n_(z+1, 0)
    REAL B_nR0(int i, int j){
      REAL dr = get_dr(i, j);
      // 2 * kappa / dr**2
      REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
      // 1 / 2 * dr
      REAL T2 = 2 * dr;
      // dk / dr
      REAL T3 = dk_dr(i, j) / dr;
      return T1 + T2 * T3;
    }

    // Calculate Coefficent for T^n_(z+1, 0)
    REAL C_nR0(int i, int j){
      REAL dr = get_dr(i, j);
      // 2 * kappa / dr**2
      REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
      // 1 / 2 * dr
      REAL T2 = 2 * dr;
      // dk / dr
      REAL T3 = dk_dr(i, j) / dr;
      return T1 - T2 * T3;
    }

    // Calculate Coefficent for T^n_(z+1, 0)
    REAL D_nR0(int i, int j){
      return D_n(i, j);
    }

    // Calculate Coefficent for T^n_(z+1, 0)
    REAL E_nR0(int i, int j){
      return E_n(i, j);
    }


  private:
    // ---------------------- PUBLIC METHODS -----------------------
    // ----------------- PRIVATE MEMBER VARIABLES  -----------------
    Field<REAL, 2> T_prime;
    int rN;
    int zN;
};
