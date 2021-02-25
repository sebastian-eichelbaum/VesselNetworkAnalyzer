//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_MATH_H
#define NOGO_MATH_H

#include <cmath>
#include <array>

#include "Types.h"

namespace nogo
{
    /**
     * Vector and 3D point representation. NOTE: std::array is a POD type -> maps to 3 consecutive floats in memory.
     */
    using Vec3 = std::array< Real, 3 >;

    /**
     * BB
     */
    struct BB2
    {
        Vec3 first;
        Vec3 second;
        //! If true, the volume is assumed to be a cylinder with z = height.
        bool isCylindric;
    };

    using BoundingBox = std::pair< Vec3, Vec3 >;

    /**
     * Calculate the length of the given vector.
     *
     * \tparam VectorType the 1x3 vector
     * \param a the vector
     *
     * \return length of the vector.
     */
    template< typename VectorType >
    typename VectorType::value_type length( const VectorType& a )
    {
        return std::sqrt( std::pow( a[ 0 ], 2 ) + std::pow( a[ 1 ], 2 ) + std::pow( a[ 2 ], 2 ) );
    }

    bool operator<=( const Vec3& a, const Vec3& b)
    {
        return ( a[ 0 ] <= b[ 0 ] ) && ( a[ 1 ] <= b[ 1 ] ) && ( a[ 2 ] <= b[ 2 ] );
    }

    /**
     * Calculate the distance between two points.
     *
     * \tparam PointType the 3 component point
     * \param a first point
     * \param b to the second point
     *
     * \return distance between both points.
     */
    template< typename PointType >
    typename PointType::value_type distance( const PointType& a, const PointType& b )
    {
        return length( PointType( { { a[ 0 ] - b[ 0 ], a[ 1 ] - b[ 1 ], a[ 2 ] - b[ 2 ] } } ) );
    }

    template< typename PointType >
    PointType abs( const PointType& a )
    {
        return PointType( { { std::abs( a[ 0 ] ), std::abs( a[ 1 ] ), std::abs( a[ 2 ] ) } } );
    }

    template< typename PointType >
    PointType difference( const PointType& a, const PointType& b )
    {
        return PointType( { { a[ 0 ] - b[ 0 ], a[ 1 ] - b[ 1 ], a[ 2 ] - b[ 2 ] } } );
    }

    template< typename PointType >
    PointType add( const PointType& a, const PointType& b )
    {
        return PointType( { { a[ 0 ] + b[ 0 ], a[ 1 ] + b[ 1 ], a[ 2 ] + b[ 2 ] } } );
    }

    //! Safely add a point b that might contain nan to another one that is safe.
    template< typename PointType >
    PointType safeAdd( const PointType& a, const PointType& b )
    {
        auto val = []( const PointType x, int i ) { return std::isnan( x[ i ] ) ? 0 : x[ i ]; };
        return PointType( { { a[ 0 ] + val( b, 0 ), a[ 1 ] + val( b, 1 ), a[ 2 ] + val( b, 2 ) } } );
    }

    template< typename PointType >
    PointType scale( const PointType& a, const Real& b )
    {
        return PointType( { { a[ 0 ] * b, a[ 1 ] * b, a[ 2 ] * b } } );
    }

    template<class T>
    constexpr const T& clamp( const T& v, const T& lo, const T& hi )
    {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }

    template< typename PointType >
    typename PointType::value_type dot( const PointType& a, const PointType& b )
    {
        return ( a[ 0 ] * b[ 0 ] ) + ( a[ 1 ] * b[ 1 ] ) + ( a[ 2 ] * b[ 2 ] );
    }

    /**
     * The normalized direction between two points from a to b.
     */
    template< typename PointType >
    PointType normalizedDirection( const PointType& a, const PointType& b )
    {
        auto diff = difference( b, a );
        return scale( diff, 1.0 / length( diff ) );
    }

    //! Calculate the angle in degree between two normalized vectors.
    template< typename PointType >
    typename PointType::value_type degree( const PointType& aNormalized, const PointType& bNormalized )
    {
        auto d = dot( aNormalized, bNormalized );
        return 180.0 * std::acos( clamp( d, -1.0, 1.0 ) ) / 3.14159265358979323846;
    }

    //! Midpoint between two given points
    template< typename PointType >
    PointType midpoint( const PointType& a, const PointType& b )
    {
        return add( a, scale( difference( b, a ), 0.5 ) );
    }
}

#endif  // NOGO_MATH_H

