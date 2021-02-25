//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_DISTANCETRANSFORM_H
#define NOGO_DISTANCETRANSFORM_H

#include <iterator>
#include <limits>

#include "Data.h"
#include "Math.h"
#include "Types.h"
#include "Range.h"

#include <vigra/multi_distance.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/multi_distance.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/multi_resize.hxx>
#include <vigra/transformimage.hxx>

#include "Logger.h"
#define LogTag "nogo/DistanceTransform"

namespace nogo
{

    /**
     * Calculate distance transform of the given volume
     *
     * \param origVolume the volume to get the transform for
     * \param calcSize the size of the volume used to do the distance transform on. 0 to use the original size.
     *
     * \return the resulting volume
     */
    SPtr< nogo::DTVolume > dt( const MaskVolume& origVolume, uint16_t calcSize = 0, uint8_t bsplineOrder = 3, uint16_t threshold = 255 )
    {
        uint16_t calcSizeX = origVolume.sizeX;
        uint16_t calcSizeY = origVolume.sizeY;
        uint16_t calcSizeZ = origVolume.sizeZ;
        float evdScale = 1.0f;

        // Prepare for scaled calculation
        if( calcSize != 0 )
        {
            // Make sure the volume is cubic
            if( ( origVolume.sizeX != origVolume.sizeY ) || ( origVolume.sizeX != origVolume.sizeZ ) )
            {
                throw std::runtime_error( "Mask volumes must have a cubic shape." );
            }

            calcSizeX = calcSize;
            calcSizeY = calcSize;
            calcSizeZ = calcSize;

            // EVD scaling factor:
            evdScale = static_cast<decltype(evdScale)>(origVolume.sizeX) / static_cast<decltype(evdScale)>(calcSizeX);
        }

        vigra::Shape3 shape( origVolume.sizeX, origVolume.sizeY, origVolume.sizeZ );
        vigra::Shape3 shapeCalc( calcSizeX, calcSizeY, calcSizeZ );
        vigra::MultiArray< 3, float > dest( shapeCalc );

        // NOTE: this is very mem-inefficient. To avoid the issue, use vigra::MultiArray in nogo::Volume.
        {
            vigra::MultiArray< 3, typename MaskVolume::value_type > source( shape );

            LogD << "Copy to source array." << LogEnd;
            for( auto z : range( origVolume.sizeZ ) )
            {
                for( auto y : range( origVolume.sizeY ) )
                {
                    for( auto x : range( origVolume.sizeX ) )
                    {
                        source( x, y, z ) = ( origVolume.data[ origVolume.index( x, y, z ) ] == 0 ) ? 0 : 255;
                    }
                }
            }

            if( shape != shapeCalc )
            {
                vigra::MultiArray< 3, typename MaskVolume::value_type > sourceCalc( shapeCalc );

                // Scale volume
                LogD << "Scaling source. B-Spline order: " << +bsplineOrder << LogEnd;
                switch (bsplineOrder)
                {
                    case 1:
                        vigra::resizeMultiArraySplineInterpolation( source, sourceCalc, vigra::BSpline<1, double>() );
                        break;
                    case 2:
                        vigra::resizeMultiArraySplineInterpolation( source, sourceCalc, vigra::BSpline<2, double>() );
                        break;
                    default:
                        vigra::resizeMultiArraySplineInterpolation( source, sourceCalc, vigra::BSpline<3, double>() );
                        break;
                }

                LogD << "Thresholding scaled source. Threshold: " << threshold << LogEnd;
                vigra::transformMultiArray( sourceCalc, sourceCalc, [threshold]( auto v ){ return  ( v > threshold ) ? 255 : 0; } );

                // Calculate Euclidean distance squared for all background pixels
                LogD << "Calculating scaled DT." << LogEnd;
                vigra::separableMultiDistance( sourceCalc, dest, true );
            }
            else
            {
                // Calculate Euclidean distance squared for all background pixels
                LogD << "Calculating DT." << LogEnd;
                vigra::separableMultiDistance( source, dest, true );
            }
        }

        LogD << "Copy to destination volume. Scaling with " << evdScale << LogEnd;
        auto dtVolume = std::make_shared< nogo::DTVolume >( calcSizeX, calcSizeY, calcSizeZ );
        for( auto z : range( calcSizeZ ) )
        {
            for( auto y : range( calcSizeY ) )
            {
                for( auto x : range( calcSizeX ) )
                {
                    dtVolume->data[ dtVolume->index( x, y, z ) ] = evdScale * dest( x, y, z );
                }
            }
        }
        return dtVolume;
    }
}

#endif // NOGO_DISTANCETRANSFORM_H


