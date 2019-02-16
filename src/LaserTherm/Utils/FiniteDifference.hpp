#ifndef LaserTherm_Utils_FiniteDifference_hpp
#define LaserTherm_Utils_FiniteDifference_hpp

/** @file FiniteDifference.hpp
  * @brief 
  * @author C.D. Clark III
  * @date 02/14/19
  */


template< typename ARRAY >
auto backDiff( const ARRAY& axis_, int i_ )
{
  if( i_ == 0 )
    return axis_[ i_ + 1 ] - axis_[ i_     ];
  else
    return axis_[ i_     ] - axis_[ i_ - 1 ];
}

template< typename ARRAY >
auto forwDiff( const ARRAY& axis_, int i_ )
{
  if(  i_ == axis_.size() - 1 )
    return axis_[ i_     ] - axis_[ i_ - 1 ];
  else
    return axis_[ i_ + 1 ] - axis_[ i_     ];
}

template< typename ARRAY >
auto centDiff( const ARRAY& axis_, int i_ )
{
  return backDiff( axis_, i_ )  + forwDiff( axis_, i_ );
}



/**
 * calculates coefficient for i_-1 term i_n central finite difference approximation of first derivative
 */
template< typename ARRAY >
auto  firstDeriv_centFDmCoeff( const ARRAY& axis_, int i_ )
{
  auto dxm = backDiff( axis_, i_ );
  auto dxc = centDiff( axis_, i_ );
  auto dxp = forwDiff( axis_, i_ );

  return -dxp / (dxc*dxm);
}

/**
 * calculates coefficient for i_ term i_n central finite difference approximation of first derivative
 */
template< typename ARRAY >
auto  firstDeriv_centFDcCoeff( const ARRAY& axis_, int i_ )
{
  auto dxm = backDiff( axis_, i_ );
  auto dxp = forwDiff( axis_, i_ );

  return (dxp - dxm) / (dxm*dxp);
}

/**
 * calculates coefficient for i_+1 term i_n central finite difference approximation of first derivative
 */
template< typename ARRAY >
auto  firstDeriv_centFDpCoeff( const ARRAY& axis_, int i_ )
{
  auto dxm = backDiff( axis_, i_ );
  auto dxc = centDiff( axis_, i_ );
  auto dxp = forwDiff( axis_, i_ );

  return dxm / (dxc*dxp);
}



/**
 * calculates coefficient for i_-1 term i_n central finite difference approximation of second derivative
 */
template< typename ARRAY >
auto  secondDeriv_centFDmCoeff( const ARRAY& axis_, int i_ )
{
  auto dxm = backDiff( axis_, i_ );
  auto dxc = centDiff( axis_, i_ );

  return 2. / (dxc*dxm);
}

/**
 * calculates coefficient for i_ term i_n central finite difference approximation of second derivative
 */
template< typename ARRAY >
auto  secondDeriv_centFDcCoeff( const ARRAY& axis_, int i_ )
{
  auto dxm = backDiff( axis_, i_ );
  auto dxp = forwDiff( axis_, i_ );

  return -2 / (dxm*dxp);
}

/**
 * calculates coefficient for i_+1 term i_n central finite difference approximation of second derivative
 */
template< typename ARRAY >
auto  secondDeriv_centFDpCoeff( const ARRAY& axis_, int i_ )
{
  auto dxp = forwDiff( axis_, i_ );
  auto dxc = centDiff( axis_, i_ );

  return 2. / (dxc*dxp);
}


#endif // include protector
