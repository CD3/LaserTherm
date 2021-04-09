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
          // r=0 case
          if(j == 0){
            //L'hospital formulas
            beta = delta_t / this->VHC[i][j];
            T1 = this->T[i][0]   * A_nBC(i , 0);
            T2 = this->T[i][1]   * B_nBC(i , 0);
            //T[i][1] from symmetry about origin
            T3 = this->T[i][1]   * C_nBC(i , 0);
            if(i == 0){
              T4 = this->minZBC.f  * D_nBC(i, 0);
              T5 = this->T[i+1][0] * E_nBC(i , 0);
            } else if(i == zN-1){
              T4 = this->T[i-1][0] * D_nBC(i , 0);
              T5 = this->maxZBC.f  * E_n(i , j);
            } else {
              T4 = this->T[i-1][0] * D_nBC(i , 0);
              T5 = this->T[i+1][0] * E_nBC(i , 0);
            }
            T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
            continue;
          }
          // r = Rmax
          if(j == rN-1){
            switch(this->maxRBC.type){
              case BC::Type::Temperature:
                beta = delta_t / this->VHC[i][j];
                T1 = this->T[i][j]   * A_n(i , j);
                T2 = this->maxRBC.f  * B_n(i , j);
                T3 = this->T[i][j-1] * C_n(i , j);
                if(i == 0){
                  T4 = this->minZBC.f  * D_n(i, j);
                  T5 = this->T[i+1][j] * E_n(i , j);
                } else if(i == zN-1){
                  T4 = this->T[i-1][j] * D_n(i , j);
                  T5 = this->maxZBC.f  * E_n(i , j);
                } else {
                  T4 = this->T[i-1][j] * D_n(i , j);
                  T5 = this->T[i+1][j] * E_n(i , j);
                }
                T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              case BC::Type::HeatFlux:
                // do HeatFlux stuff
                break;
              default:
                throw 42;
            }
            continue;
          }
          // z = 0
          if(i == 0){
            switch(this->minZBC.type){
              case BC::Type::Temperature:
                beta = delta_t / this->VHC[i][j];
                T1 = this->T[0][j]   * A_n(i , j);
                T2 = this->T[0][j+1] * B_n(i , j);
                T3 = this->T[0][j-1] * C_n(i , j);
                T4 = this->minZBC.f  * D_n(i , j);
                T5 = this->T[1][j] * E_n(i , j);
                T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              case BC::Type::HeatFlux:
                // do HeatFlux stuff
                break;
              default:
                throw 42;
            }
            continue;
          }
          // z = Zmax
          if(i == zN-1){
            switch(this->maxZBC.type){
              case BC::Type::Temperature:
                beta = delta_t / this->VHC[i][j];
                T1 = this->T[i][j]   * A_n(i , j);
                T2 = this->T[i][j+1] * B_n(i , j);
                T3 = this->T[i][j-1] * C_n(i , j);
                T4 = this->T[i-1][j] * D_n(i , j);
                T5 = this->maxZBC.f  * E_n(i , j);
                T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
                break;
              case BC::Type::None:
                // do stuff for temp type
                break;
              case BC::Type::HeatFlux:
                // do HeatFlux stuff
                break;
              default:
                throw 42;
            }
            continue;
          }
          beta = delta_t / this->VHC[i][j];
          T1 = this->T[i][j]   * A_n(i , j);
          T2 = this->T[i][j+1] * B_n(i , j);
          T3 = this->T[i][j-1] * C_n(i , j);
          T4 = this->T[i-1][j] * D_n(i , j);
          T5 = this->T[i+1][j] * E_n(i , j);
          T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
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
    // Finite Difference deriv for kappa wrt r
    REAL dk_dr(int i, int j){
      REAL numer = (this->k[i][j + 1] - this->k[i][j - 1]);
      // is this an ok way to do 2\Delta r?
      REAL denom = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j - 1)[1];
      return numer / denom;
    }
    // Finite Difference deriv for kappa wrt z
    REAL dk_dz(int i, int j){
      REAL numer = (this->k[i + 1][j] - this->k[i - 1][j]);
      // is this an ok way to do 2\Delta r?
      REAL denom = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i - 1, j)[0];
      return numer / denom;
    }
    // Calculate Coefficent for T^n_(r, z)
    REAL A_n(int i, int j){
      // k.getCoord returns an array, hopefully like {z, r} to be consistent
      // so (r, z) is written as (k.getCoord[1], k.getCoord[0])
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // -2 * kappa / dr**2
      REAL T1 = -2 * this->k[i][j] / pow(dr, 2);
      // -2 * kappa / dz**2
      REAL T2 = -2 * this->k[i][j] / pow(dz, 2);
      return T1 + T2;
    }
    REAL A_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      REAL T1 = (-4 * this->k[i][j]) / (dr * dr);
      REAL T2 = (-2 * this->k[i][j]) / (dz * dz);
      return T1 + T2;
    }
    // Calculate Coefficent for T^n_(r+1, z)
    REAL B_n(int i, int j){
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // kappa / 2 r dr
      REAL T1 = this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      REAL T2 = this->k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      REAL T3 = dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }
    REAL B_nBC(int i, int j){
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      // 2 * kappa / dr**2
      REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
      // 1 / 2 * dr
      REAL T2 = 2 * dr;
      // dk / dr
      REAL T3 = dk_dr(i, j) / dr;
      return T1 + T2 * T3;
    }
    // Calculate Coefficent for T^n_(r-1, z)
    REAL C_n(int i, int j){
      //return (N) 0;
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // - kappa / 2r * dr
      REAL T1 = - this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      REAL T2 = this->k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      REAL T3 = - dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }
    REAL C_nBC(int i, int j){
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      // 2 * kappa / dr**2
      REAL T1 = 2 * this->k[i][j] / pow(dr, 2);
      // 1 / 2 * dr
      REAL T2 = 2 * dr;
      // dk / dr
      REAL T3 = dk_dr(i, j) / dr;
      return T1 - T2 * T3;
    }
    // Calculate Coefficent for T^n_(r, z-1)
    REAL D_n(int i, int j){
      //return (N) 0;
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // kappa / dz**2
      REAL T1 = this->k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      REAL T2 = - dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }
    REAL D_nBC(int i, int j){
      return D_n(i, j);
    }
    // Calculate Coefficent for T^n_(r, z+1)
    REAL E_n(int i, int j){
      //return (N) 0;
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // kappa / dz**2
      REAL T1 = this->k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      REAL T2 = dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }
    REAL E_nBC(int i, int j){
      return E_n(i, j);
    }


  private:
    // ---------------------- PUBLIC METHODS -----------------------
    // ----------------- PRIVATE MEMBER VARIABLES  -----------------
    Field<REAL, 2> T_prime;
    int rN;
    int zN;
};
