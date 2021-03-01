#include <libField/Field.hpp>
#include "LaserTherm/HeatSolvers/BoundaryConditions.hpp"

namespace BC = HeatSolvers::BoundaryConditions;

template <class N>
class Explicit{
  public:
    // -------------------- PUBLIC CONSTRUCTORS --------------------
    Explicit(int _rj, int _zi){
      this->rj = _rj;
      this->zi = _zi;
      T.reset(zi, rj);
      A.reset(zi, rj);
      VHC.reset(zi, rj);
      k.reset(zi, rj);
    }
    Explicit(){
      // Default constructor to be safe
    }
    // ---------------------- PUBLIC METHODS -----------------------
    // delta_t is a step forward in time in seconds
    void stepForward(N delta_t){
      N beta;
      Field<N, 2> T_prime(zi, rj);
      for(N t = 0.0; t < delta_t; t+=dt){
        for(int i = 0; i < zi; i++){
          for(int j = 0; j < rj; j++){
            beta = delta_t / VHC[i][j];
            T_prime[i][j] = T[i][j] + beta
            * (T[i][j] * A_n(i , j)
            + T[i][j+1] * B_n(i , j)
            + T[i][j-1] * C_n(i , j)
            + T[i-1][j] * D_n(i , j)
            + T[i+1][j] * E_n(i , j));
          }
        }
        T.set(1.0);
        T *= T_prime;
      }
    }
    int get_rdims(){
      return this->rj;
    }
    int get_zdims(){
      return this->zi;
    }
    // ----------------- PUBLIC MEMBER VARIABLES -------------------
    Field<N, 2> T;
    Field<N, 2> A;
    Field<N, 2> VHC;
    Field<N, 2> k;
    double dt = 0.001;

  /* Protected namespace is allowed to be
   * accessed by classes that inhereit using
   * public or protected keywords, but not by
   * outside or unrelated classes */
  protected:
    // ---------------- PROTECTED MEMBER VARIABLES  ----------------
    // --------------------- PROTECTED METHODS ---------------------
    // Calculate Coefficent for T^n_(r, z)
    N A_n(int i, int j){

    }
    N A_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
      // double fprime = BC::FiniteDifference<N>::f;
    }
    // Calculate Coefficent for T^n_(r+1, z)
    N B_n(int i, int j){

    }
    N B_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r-1, z)
    N C_n(int i, int j){

    }
    N C_nBC(int i, int j){
    // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r, z-1)
    N D_n(int i, int j){

    }
    N D_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }
    // Calculate Coefficent for T^n_(r, z+1)
    N E_n(int i, int j){

    }
    N E_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
    }


  private:
    // ---------------------- PUBLIC METHODS -----------------------
    // ----------------- PRIVATE MEMBER VARIABLES  -----------------
    int rj;
    int zi;
};
