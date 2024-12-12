#pragma once

/** @file Slab.hpp
 * @brief
 * @author C.D. Clark III
 * @date 03/12/19
 */

#include <optional>

namespace Structures::_1D::Cartesian
{
template<typename REAL>
class Slab
{
 public:
  std::optional<REAL> minSurfacePosition;
  std::optional<REAL> maxSurfacePosition;

  Slab() = default;

  void setMinSurfacePosition(REAL p) { minSurfacePosition = p; }

  void setMaxSurfacePosition(REAL p) { maxSurfacePosition = p; }

  REAL getMinSurfacePosition() const { return minSurfacePosition.value(); }
  REAL getMaxSurfacePosition() const { return maxSurfacePosition.value(); }
  REAL getThickness() const
  {
    return maxSurfacePosition.value() - minSurfacePosition.value();
  }

  void setThickness(REAL t)
  {
    if (minSurfacePosition) maxSurfacePosition = minSurfacePosition.value() + t;
    if (maxSurfacePosition) minSurfacePosition = maxSurfacePosition.value() - t;
  }

  bool isInside(REAL x)
  {
    bool inside = false;

    if (minSurfacePosition && maxSurfacePosition) {
      if (minSurfacePosition <= x && x <= maxSurfacePosition) inside = true;
    }
    if (!minSurfacePosition && maxSurfacePosition) {
      if (x <= maxSurfacePosition) inside = true;
    }
    if (minSurfacePosition && !maxSurfacePosition) {
      if (minSurfacePosition <= x) inside = true;
    }

    return inside;
  }
};

}  // namespace Structures::_1D::Cartesian
