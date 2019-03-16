#ifndef LaserTherm_Structures__1D_Slab_hpp
#define LaserTherm_Structures__1D_Slab_hpp

/** @file Slab.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/12/19
  */


namespace Structure::_1D
{
template<class MaterialType>
class Slab
{
  public:
    MaterialType material;

    Slab( const MaterialType mat ):material(mat){}
    Slab() = default;

    void setMaterial(const MaterialType mat )
    {
      material = mat;
    }

};

}


#endif // include protector
