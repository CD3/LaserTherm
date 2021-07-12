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
      this->T_prime.set(0);
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
    
    void stepForward(REAL delta_t, int n){
      for(int i = 0; i < n i++){
        this->stepForward(delta_t);
      }
    }
    
    /** T0 - initial temperature solution, T_Dt - the temperature solution at a
     * time Dt in the future obtained by a single application of the heat
     * solver, T_Dt2 - the temperature solution at a time Dt/2 in the future,
     * T2_Dt - the temperature solution at a time Dt in the future
     * obtained by two applications of the heat solver based on half time steps
     *
     * time step optimization occurs here 
     * 
     * T_Dt is determined WITHIN the call,
     * so auto convergence flag must already be called, possible change later 
     *
     * The algorithm is kind of about BOTH calculating the minimum convergent
     * timestep AND moving forward by by two steps of that timestep and that
     * kind of just feels weird... there should be a seperate functionality to
     * FIND the minimum time step and a seperate function to actually move
     * forward.
     * Currently I'm thinking
     * - function to find min convergent time step
     * - overloaded stepforward(REAL a, REAL b) @param a the total time to move
     *   by @param b the increment (step size) to use (assert a > b) 
     * - function to JUST DO IT (moveForward a macro Deltat without returning
     *   the timestep used), results in a loss of information but may be
     *   valuable for a simpler move
     *
     * */
    REAL findTimeStep(REAL Delta_t, REAL t_min){
      // use T0 store initial state, so this-> MUST be initialized with temperature at t=0
      bool invalidNodeFound = false;
      Field<REAL, 2> T0, T_Dt, , T_Dt2, T2_Dt;
      T0.reset(this->T.getCoordinateSystemPtr());
      T_Dt.reset(this->T.getCoordinateSystemPtr());
      T_Dt2.reset(this->T.getCoordinateSystemPtr());
      T2_Dt.reset(this->T.getCoordinateSystemPtr());

      T0.set(0);
      // T0 is holding T0, the initial temperature solution
      T0 += this->T;
      // Calculate T_Dt (only has to happen once)
      this->stepForward(Delta_t);
      T_Dt.set(0);
      T_Dt += this->T;
      while(invalidNodeFound){
        T_Dt2.set(0);
        T2_Dt.set(0);
        // reset T to original value to prep for next calculation
        this->T.set(0);
        this->T += T0;
        // calculate T2_dt, and store T_dt2
        this->stepForward(Delta_t/2);
        T_Dt2 += this->T;
        this->stepForward(Delta_t/2);
        T2_Dt += this->T;
        //double totalErr(Field<REAL, 2>& f, REAL (*sol)(REAL, REAL, REAL), REAL (*err)(REAL, REAL), REAL t){
        for(int i = 0; i < this->T.size(0); i++){
          for(int j = 0; j < this->T.size(1); j++){
            invalidNodeFound = totalErr() > 100;
            invalidNodeFound ? break break : ;
            if(invalidNodeFound){
              // Replace T_Dt with T_Dt2 ???
              T_Dt.set(0);
              T_Dt += T_Dt2;
              if(Delta_t / 2 <= t_min){
                return t_min;
              }
              Delta_t /= 2;
              break;
            }
          }
          if(invalidNodeFound) { break; }
          // LASTMARK
          // find some way to code 3(b)(i) with conditionals...
        }
      }
      // calculate single initial guess for timestep, save all or some of the values
      // calculate 2 t/2 timesteps
      // compare error
      // rinse n repeat
      // once have smallest time step, t*
      // divide Delta_t by t*, and run
      // stepForward(t*, Delta_t / t*)
      return Delta_t / 2;
      // lol
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
    double totalErr(){
      double e = 0.0;
    }
};
