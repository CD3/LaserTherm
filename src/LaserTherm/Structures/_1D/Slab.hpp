#ifndef LaserTherm_Structures__1D_Slab_hpp
#define LaserTherm_Structures__1D_Slab_hpp

/** @file Slab.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 03/12/19
  */

#include <optional>

namespace Structures::_1D
{
template<class MaterialType, class REAL>
class Slab
{
  public:

    std::optional<REAL> minSurfacePosition;
    std::optional<REAL> maxSurfacePosition;
    std::optional<MaterialType> material;

    Slab( const MaterialType mat ):material(mat){}
    Slab() = default;

    void setMaterial(const MaterialType mat )
    {
      material = mat;
    }

    void setMinSurfacePosition(REAL p)
    {
      minSurfacePosition = p;
    }

    void setMaxSurfacePosition(REAL p)
    {
      maxSurfacePosition = p;
    }

    void setThickness(REAL t)
    {
      if( minSurfacePosition )
        maxSurfacePosition = minSurfacePosition.get() + t;
      if( maxSurfacePosition )
        minSurfacePosition = maxSurfacePosition.get() - t;
    }

    bool isInside( REAL x )
    {
      bool inside = false;

      if( minSurfacePosition && maxSurfacePosition )
      {
        if( minSurfacePosition <= x && x <= maxSurfacePosition )
          inside = true;
      }
      if( !minSurfacePosition && maxSurfacePosition )
      {
        if( x <= maxSurfacePosition )
          inside = true;
      }
      if( minSurfacePosition && !maxSurfacePosition )
      {
        if( minSurfacePosition <= x )
          inside = true;
      }

      return inside;
    }

};

}


#endif // include protector
