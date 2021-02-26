#include <libField/Field.hpp>

template <class N>
class Explicit{
  public:
    // constructors
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
    // methods
    void stepForward(N delta_t){
      // junk goes here wow i spelled everythign right
    }
    // Calculate Coefficent for T^n_(r, z)
    N A_n(int i, int j){

    }
    N A_nBC(int i, int j){
      // Some reference to Boundary Conditions set elsewhere
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
    // member variables
    int rj;
    int zi;
    Field<N, 2> T;
    Field<N, 2> A;
    Field<N, 2> VHC;
    Field<N, 2> k;
  protected:
  private:
};
