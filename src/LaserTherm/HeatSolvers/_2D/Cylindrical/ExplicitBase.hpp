#include <libField/Field.hpp>
#include <math.h>
#include <exception>
#include <stdexcept>
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
          zMax = i == zN-1; zMin = i == 0; // coefficient for calculations
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
                throw std::runtime_error(std::string("No BC Type!"));
            }
          } else {
            // if not on r boundary give 'normal T2 T3'
            T2 = this->T[i][j-1] * static_cast<IMP*>(this)->B_n(i , j);
            T3 = this->T[i][j+1] * static_cast<IMP*>(this)->C_n(i , j); } if(zMin){
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
                throw std::runtime_error(std::string("No BC Type!"));
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
                throw std::runtime_error(std::string("No BC Type!"));
            }
          } else {
            // if not on z boundary give 'normal T4 T5'
            T4 = this->T[i-1][j] * static_cast<IMP*>(this)->D_n(i , j);
            T5 = this->T[i+1][j] * static_cast<IMP*>(this)->E_n(i , j);
          }
          this->T_prime[i][j] = beta * (T1 + T2 + T3 + T4 + T5 + this->A[i][j]);
          if(isnan(this->T_prime[i][j])){
            // make this a better error
            throw std::runtime_error(std::string("NaN Encountered!"));
          }
        }
      }
      this->T += T_prime;
    }
    
    /**
     * a function to evolve this->T in time a total of delta_t, in increments of dt
     *
     * @param delta_t the desired time to 'move forward' in one step
     * @param n the number of steps to take
     * */
    void stepForward(REAL delta_t, int n){
      for(int i = 0; i < n; i++){
        this->stepForward(delta_t);
      }
    }
    
    /**
     * a function to evolve this->T in time a total of delta_t, in increments of dt
     *
     * @param delta_t the total desired time to 'move forward' (non inclusive)
     * @param dt the increment taken to reach delta_t
     * */
    void stepForward(REAL delta_t, REAL dt){
      assert(delta_t > dt);
      this->stepForward(dt, static_cast<int>(delta_t / dt));
    }

    /** 
     * A method (roughly) implementing the Automatic Converge Time Step algorithm as detailed in /doc/writups/...
     *
     * @param Delta_t is the 'initial guess' for the minimum convergent time step
     * @param t_min is the minimum time step allowed by the algorithm, 0 by default
     *
     * @returns A REAL representing an optimized timestep where the evolution of
     * T is convergent
     *
     * */
    REAL findTimeStep(REAL Delta_t, REAL t_min=0.0){
      assert(t_min >= 0); 
      assert(Delta_t > 0); 

      bool invalidNodeFound = true;
      //T0 = T(0), T_Dt = T(\Delta t), T_Dt2 = T(\Delta t / 2), T2_Dt = T^2(\Delta t)
      Field<REAL, 2> T0, T_Dt, T_Dt2, T2_Dt;

      T0.reset(this->T.getCoordinateSystemPtr());
      T_Dt.reset(this->T.getCoordinateSystemPtr());
      T_Dt2.reset(this->T.getCoordinateSystemPtr());
      T2_Dt.reset(this->T.getCoordinateSystemPtr());

      // This set(0), += method is a shorthand I found to 'copy' a Field's data
      // to another without worrying about overwriting data or passing a
      // reference by accident, and I use it a lot
      // Make T0 have the same profile as T
      T0.set(0);
      T0 += this->T;
      // Calculate and set T_Dt (calculation only has to happen once)
      this->stepForward(Delta_t);
      T_Dt.set(0);
      T_Dt += this->T;

      // Enter algorithm
      while(invalidNodeFound){

        // At the beginning of each iteration, clear values
        T_Dt2.set(0);
        T2_Dt.set(0);
        // reset T to initial value to prep for next calculation since it
        // changes with each stepForward
        this->T.set(0);
        this->T += T0;
        // calculate and store T_Dt2, and store T_Dt2
        this->stepForward(Delta_t/2);
        T_Dt2 += this->T;
        // calculate and store T2_Dt, and store T2_Dt
        this->stepForward(Delta_t/2);
        T2_Dt += this->T;

        // I'm parsing through the field here but only really using the 'average method'...
        for(int i = 0; i < T0.size(0); i++){
          for(int j = 0; j < T0.size(1); j++){
            // set to true if error limit is surpassed, false otherwise
            invalidNodeFound = randSampleErr(T_Dt, T2_Dt, percErr, 10) > 0.01;
            // invalidNodeFound = avgErr(T_Dt, T2_Dt, percErr) > 0.001;
            // invalidNodeFound = totalErr(T_Dt, T2_Dt, percErr) > 1;
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
        }
      }
      this->T.set(0);
      this->T += T0;
      return Delta_t;
    }


    /**
     * A function to evolve this->T in time a total of Delta_t, in increments
     * based on an optimized time step larger than some minimum
     *
     * @param Delta_t the desired time to 'move forward' in one step 
     * @param t_min the minimum time step allowed to increment by (0 by default)
     * */
    REAL moveForward(REAL Delta_t, REAL t_min=0.0){
      assert(Delta_t > t_min);
      REAL dt = this->findTimeStep(Delta_t, t_min);
      stepForward(Delta_t, dt); 
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
    static double percErr(REAL num, REAL correct){
      return abs((num - correct) / (num + correct));
    }

    double randSampleErr(Field<REAL, 2>& T1, Field<REAL, 2>& T2, double (*err)(REAL, REAL), int n){
      assert(T1.size(0) == T2.size(0));
      assert(T1.size(1) == T2.size(1));
      double e = 0;
      int zi, rj;
      int rMAX = T1.size(1);
      int zMAX = T1.size(0);
      for(int i = 0; i < n; i++){
        rj = rand() % rMAX;
        zi = rand() % zMAX;
        e += err(T1[zi][rj], T2[zi][rj]);
      }
      return e / n;
    }

    double avgErr(Field<REAL, 2>& T1, Field<REAL, 2>& T2, double (*err)(double, double)){
      assert(T1.size(0) == T2.size(0));
      assert(T1.size(1) == T2.size(1));
      double e = 0;
      int rMAX = T1.size(1);
      int zMAX = T1.size(0);
      for(int i = 0; i < zMAX; i++){
        for(int j = 0; j < rMAX; j++){
          e += err(T1[i][j], T2[i][j]);
        }
      }
      e /= zMAX * rMAX;
      return e;
    }

    double totalErr(Field<REAL, 2>& T1, Field<REAL, 2>& T2, double (*err)(double, double)){
      assert(T1.size(0) == T2.size(0));
      assert(T1.size(1) == T2.size(1));
      double e = 0;
      int rMAX = T1.size(1);
      int zMAX = T1.size(0);
      for(int i = 0; i < zMAX; i++){
        for(int j = 0; j < rMAX; j++){
          e += err(T1[i][j], T2[i][j]);
        }
      }
      return e;
    }
};
