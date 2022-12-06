//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#include <iostream>

#include "Data.h"
#include "IO.h"
#include "Analyze.h"
#include "DistanceTransform.h"

#include <vigra/multi_array.hxx>
#include <vigra/resizeimage.hxx>
#include <vigra/multi_resize.hxx>

#include "Logger.h"
#define LogTag "prepareCT"

namespace
{
    /**
     * Print welcome message and copyright information.
     */
    void printHello()
    {
        std::cout << "NOGO Prepare CT - part of the NOGO Result Analyzer Software Package" << std::endl
                  << "Copyright (C) 2015  Sebastian Eichelbaum (http://www.nemtics.com)" << std::endl
                  << std::endl
                  << "This program comes with ABSOLUTELY NO WARRANTY." << std::endl
                  << "For more details, refer to the license that" << std::endl
                  << "has been shipped along with this program." << std::endl
                  << std::endl;
    }

    /**
     * Print usage info.
     */
    void printUsage()
    {
        std::cout << "Usage:" << std::endl << std::endl
                  << "preparect OUTPUT_PATH VOLUME_FILE" << std::endl << std::endl
                  << "OUTPUT_PATH - the relative or absolute path where to write results. Needs to exist." << std::endl
                  << "VOLUME_FILE - the relative or absolute path to the volume file." << std::endl;
    }
}

/**
 * The program's main function and entry point. Start reading here.
 *
 * \param argc number of arguments given on command line
 * \param argv argument list
 *
 * \return 0 if everything was fine. Return something else to denote an error.
 */
int main( int argc, char** argv )
{
    printHello();

    // Parameter validity
    if( argc != 3 )
    {
        printUsage();
        return -1;
    }

    LogI << "Startup successful. Starting ..." << LogEnd;
    LogD << "Parameter Summary:" << LogEnd;
    for( auto argumentIndex : nogo::range( 0, argc ) )
    {
        LogD << argumentIndex << ": " << argv[ argumentIndex ] << LogEnd;
    }

    // Catch any kind of exception
    try
    {
        // 1 - Loading the data
        auto volumeFile = std::string( argv[ 2 ] );
        auto outputDir = std::string( argv[ 1 ] );

        if( outputDir.empty() )
        {
            throw std::runtime_error( "Output path is not specified." );
        }
        if( ( outputDir.back() != '/' ) &&  ( outputDir.back() != '\\' ) )
        {
            outputDir += "/";
        }

        LogI << "Try loading volume from \"" << volumeFile << "\"." << LogEnd;
        auto ctVolume = std::make_shared< nogo::CTVolume >( 512, 512, 512 );
        nogo::loadVolume( volumeFile, ctVolume );
        auto shrunkVolume = std::make_shared< nogo::Volume< float > >( 512, 512, 512 );

        // NOTE: this is very mem-inefficient. To avoid the issue, use vigra::MultiArray in nogo::Volume!? Involves work. No time yet.
        {
            vigra::Shape3 shape( ctVolume->sizeX, ctVolume->sizeY, ctVolume->sizeZ );
            vigra::Shape3 shrunkShape( 512, 512, 512 );

            vigra::MultiArray< 3, typename nogo::CTVolume::value_type > source( shape );
            vigra::MultiArray< 3, float > shrunk( shrunkShape );

            LogD << "Copy to source array." << LogEnd;
            for( auto z : nogo::range( ctVolume->sizeZ ) )
            {
                for( auto y : nogo::range( ctVolume->sizeY ) )
                {
                    for( auto x : nogo::range( ctVolume->sizeX ) )
                    {
                        // mask
                        auto p = nogo::Vec3( { { static_cast< nogo::Real >( x ), static_cast< nogo::Real >( y ), static_cast< nogo::Real >( z ) } } );
                        int isIn = nogo::pointInVolume( p, nogo::cutBoundingBox( ctVolume->boundingBox, 98, 100.0 ), true ) ? 1 : 0;

                        // Write
                        source( x, y, z ) = isIn ? ctVolume->data[ ctVolume->index( x, y, z ) ] : 0;
                    }
                }
            }

            // Scale
            LogD << "Scale." << LogEnd;
            vigra::resizeMultiArraySplineInterpolation( source, shrunk, vigra:: CatmullRomSpline< double >() );

            LogD << "Copy to destination volume." << LogEnd;

            // Copy to target. Reuse original mem
            for( auto z : nogo::range( ctVolume->sizeZ ) )
            {
                for( auto y : nogo::range( ctVolume->sizeY ) )
                {
                    for( auto x : nogo::range( ctVolume->sizeX ) )
                    {
                        ctVolume->data[ ctVolume->index( x, y, z ) ] = source( x, y, z );
                    }
                }
            }

            // Copy to target. Reuse original mem
            for( auto z : nogo::range( shrunkVolume->sizeZ ) )
            {
                for( auto y : nogo::range( shrunkVolume->sizeY ) )
                {
                    for( auto x : nogo::range( shrunkVolume->sizeX ) )
                    {
                        shrunkVolume->data[ shrunkVolume->index( x, y, z ) ] = shrunk( x, y, z );
                    }
                }
            }

        }
        nogo::write( outputDir + "/ctMasked.raw", *ctVolume );
        nogo::write( outputDir + "/ctMaskedShrunk.raw", *shrunkVolume );

        LogI << "Analysis done. Quit." << LogEnd;
    }
    catch( std::exception& e )
    {
        LogE << "Exception caught: " << e.what() << LogEnd;
        return -1;
    }

    LogI << "Bye." << LogEnd;

    return 0;
}


