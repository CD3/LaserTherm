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
    void A_n(int i, int j){

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
