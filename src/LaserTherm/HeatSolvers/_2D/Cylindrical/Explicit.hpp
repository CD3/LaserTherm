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
    }
    Explicit(){
      // Default constructor to be safe
    }
    // ---------------------- PUBLIC METHODS -----------------------
    // delta_t is a step forward in time in seconds
    void stepForward(REAL delta_t){
      REAL beta;
      Field<REAL, 2> T_prime(zN, rN);
      for(int i = 1; i < zN - 1; i++){
        for(int j = 1; j < rN - 1; j++){
          beta = delta_t / this->VHC[i][j];
          REAL T1 = this->T[i][j]   * A_n(i , j);
          REAL T2 = this->T[i][j+1] * B_n(i , j);
          REAL T3 = this->T[i][j-1] * C_n(i , j);
          REAL T4 = this->T[i-1][j] * D_n(i , j);
          REAL T5 = this->T[i+1][j] * E_n(i , j);
          T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5);
          if(isnan(T_prime[i][j])){
            // make this a better error
            throw 20;
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
      // double fprime = BC::FiniteDifference<REAL>::f;
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
      // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r-1, z)
    REAL C_n(int i, int j){
      //return (N) 0;
      REAL dr = this->k.getCoord(i, j + 1)[1] - this->k.getCoord(i, j)[1];
      REAL dz = this->k.getCoord(i + 1, j)[0] - this->k.getCoord(i, j)[0];
      // -this->kappa / 2r * dr
      REAL T1 = - this->k[i][j] / (2 * this->k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      REAL T2 = this->k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      REAL T3 = - dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }
    REAL C_nBC(int i, int j){
    // Some reference to Boundary Conditions set elsewhere
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
      // Some reference to Boundary Conditions set elsewhere
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
      // Some reference to Boundary Conditions set elsewhere
    }


  private:
    // ---------------------- PUBLIC METHODS -----------------------
    // ----------------- PRIVATE MEMBER VARIABLES  -----------------
    int rN;
    int zN;
};
