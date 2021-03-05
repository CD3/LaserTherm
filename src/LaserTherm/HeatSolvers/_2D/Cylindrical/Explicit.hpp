#include <libField/Field.hpp>
#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"

namespace BC = HeatSolvers::BoundaryConditions;

template <class N>
class Explicit{
  public:
    // -------------------- PUBLIC CONSTRUCTORS --------------------
    Explicit(int _rN, int _zN){
      this->rN = _rN;
      this->zN = _zN;
      T.reset(zN, rN);
      A.reset(zN, rN);
      VHC.reset(zN, rN);
      k.reset(zN, rN);
    }
    Explicit(){
      // Default constructor to be safe
    }
    // ---------------------- PUBLIC METHODS -----------------------
    // delta_t is a step forward in time in seconds
    void stepForward(N delta_t){
      N beta;
      Field<N, 2> T_prime(zN, rN);
      for(int i = 1; i < zN - 1; i++){
        for(int j = 1; j < rN - 1; j++){
          beta = delta_t / VHC[i][j];
          // probably make this a function? vvv
          // Sink boundary conditions (if we go from 1 to )
          // Make Tprime the increase amount and just add it
          T_prime[i][j] = beta
          * (T[i][j] * A_n(i , j)
          + T[i][j+1] * B_n(i , j)
          + T[i][j-1] * C_n(i , j)
          + T[i-1][j] * D_n(i , j)
          + T[i+1][j] * E_n(i , j));
        }
      }
      T += T_prime;
    }
    int get_rdims(){
      return this->rN;
    }
    int get_zdims(){
      return this->zN;
    }
    // ----------------- PUBLIC MEMBER VARIABLES -------------------
    Field<N, 2> T;
    Field<N, 2> A;
    Field<N, 2> VHC;
    Field<N, 2> k;
    // we (probably) can't explicitly define dt as N type :(
    double dt = 0.001;

  protected:
    // ---------------- PROTECTED MEMBER VARIABLES  ----------------
    // --------------------- PROTECTED METHODS ---------------------
    // Finite Difference deriv for kappa wrt r
    N dk_dr(int i, int j){
      N numer = (k[i][j + 1] - k[i][j - 1]);
      // is this an ok way to do 2\Delta r?
      N denom = k.getCoord(i, j + 1)[1] - k.getCoord(i, j - 1)[1];
      return numer / denom;
    }
    // Finite Difference deriv for kappa wrt z
    N dk_dz(int i, int j){
      N numer = (k[i + 1][j] - k[i - 1][j]);
      // is this an ok way to do 2\Delta r?
      N denom = k.getCoord(i + 1, j)[0] - k.getCoord(i - 1, j)[0];
      return numer / denom;
    }
    // Calculate Coefficent for T^n_(r, z)
    N A_n(int i, int j){
      // k.getCoord returns an array, hopefully like {z, r} to be consistent
      // so (r, z) is written as (k.getCoord[1], k.getCoord[0])
      N dr = k.getCoord(i, j + 1)[1] - k.getCoord(i, j)[1];
      N dz = k.getCoord(i + 1, j)[0] - k.getCoord(i, j)[0];
      // -2 * kappa / dr**2
      N T1 = -2 * k[i][j] / pow(dr, 2);
      // -2 * kappa / dz**2
      N T2 = -2 * k[i][j] / pow(dz, 2);
      return T1 + T2;
    }
    N A_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
      // double fprime = BC::FiniteDifference<N>::f;
    }
    // Calculate Coefficent for T^n_(r+1, z)
    N B_n(int i, int j){
      N dr = k.getCoord(i, j + 1)[1] - k.getCoord(i, j)[1];
      N dz = k.getCoord(i + 1, j)[0] - k.getCoord(i, j)[0];
      // kappa / 2 r dr
      N T1 = k[i][j] / (2 * k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      N T2 = k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      N T3 = dk_dr(i, j) / (2 * dr);
      return T1 + T2 + T3;
    }
    N B_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r-1, z)
    N C_n(int i, int j){
      //return (N) 0;
      N dr = k.getCoord(i, j + 1)[1] - k.getCoord(i, j)[1];
      N dz = k.getCoord(i + 1, j)[0] - k.getCoord(i, j)[0];
      // -kappa / 2r * dr
      N T1 = - k[i][j] / (2 * k.getCoord(i, j)[1] * dr);
      // kappa / dr**2
      N T2 = k[i][j] / pow(dr, 2);
      // dk_dr / 2 dr
      N T3 = - dk_dr(i, j) / (2 * dr); 
      return T1 + T2 + T3;
    }
    N C_nBC(int i, int j){
    // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r, z-1)
    N D_n(int i, int j){
      //return (N) 0;
      N dr = k.getCoord(i, j + 1)[1] - k.getCoord(i, j)[1];
      N dz = k.getCoord(i + 1, j)[0] - k.getCoord(i, j)[0];
      // kappa / dz**2
      N T1 = k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      N T2 = - dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }
    N D_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r, z+1)
    N E_n(int i, int j){
      //return (N) 0;
      N dr = k.getCoord(i, j + 1)[1] - k.getCoord(i, j)[1];
      N dz = k.getCoord(i + 1, j)[0] - k.getCoord(i, j)[0];
      // kappa / dz**2
      N T1 = k[i][j] / pow(dz, 2);
      // - dz_dr / 2 dz
      N T2 = dk_dz(i, j) / (2 * dz);
      return T1 + T2;
    }
    N E_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }


  private:
    // ---------------------- PUBLIC METHODS -----------------------
    // ----------------- PRIVATE MEMBER VARIABLES  -----------------
    int rN;
    int zN;
};
