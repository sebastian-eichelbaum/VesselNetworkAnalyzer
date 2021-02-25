//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#pragma once

#include <string>
#include <fstream>
#include <streambuf>
#include <exception>

#include <picojson.h>

#include "Types.h"

#include "Logger.h"
#define LogTag "nogo/Settings"

namespace nogo
{
    struct Settings
    {
        // The edge size of the volume. Raw volumes have to have this size cubed.
        int volumeSize = 1024;

        // Scale used for the EVD. The EVD is calculated in voxel space. Use the voxel size in the volume as scaler.
        Real pixelSize = 0.73;

        // If true, the network data is in voxel space and will be scaled by pixel size
        bool networkIsVoxelSpace = false;

        // If true, the volume will be seen as a cylinder instead of a cube.
        bool cylindric = true;

        // If true, the whole network is assumed to be cropped already. It will be used as is but the volume will be
        // assumed to be defined by networkCutXY and networkCutZ (i.e. for densities).
        bool preCroppedNetwork = false;

        // Data can be cut using a cylinder. This defines the diameter of the cylinder in percent of the original
        // volume and its height in percent of the original volume.
        Real networkCutXY = 50.0;
        Real networkCutZ = 75.0;

        // Like the network data, the masked volume can be cut in a similar way.
        Real volumeCutXYVVF = 95.0;
        Real volumeCutZVVF = 100.0;

        // Like the network data, the masked volume can be cut in a similar way.
        Real volumeCutXYEVD = 100.0;
        Real volumeCutZEVD = 100.0;

        // Restrict the max degree of branch points. Useful to filter some noisy branch points
        bool restrictDegree = true;
        // If restricted, map larger branch point degree to the max?
        bool mapLargerToMax = false;
        // If restricted, use this a s max branch point degree.
        int maxBranchpointDegree = 6;

        // Should very close branch points be merged into one?
        bool mergeCloseBranchPoints = true;
        // If merging close branch points, use this as a max distance to merge in the same coordinate system as the
        // network data.
        Real maxMergeDistance = 2.0;

        // We store the diameter per segment length in a discrete length-span histogram.
        Real diameterPerLengthMin = 1;
        Real diameterPerLengthMax = 70;
        int diameterPerLengthNumBins = 14;

        // Force a min and max range for EVD values.
        Real minEVD = 0.10;
        bool useMinEVD = true;
        Real maxEVD = 200.0;
        bool useMaxEVD = true;

        // Define a discretization scheme for EVD. This defines the histogram used to categorize EVD values.
        Real evdHistMin = 0;
        Real evdHistMax = 200.0;
        int evdNumBins = 20;

        // Scale the volume for evd calculation and saving. Set to 0 (default if unset) to disable any scaling. When
        // changing this value, old saved EVD volumes are probably not useable anymore.
        int evdVolumeSize = 0;

        // When scaling the EVD volume, use this threshold to define the vessel mask threshold in the scaled volume.
        int evdScaleThreshold = 192;

        // When scaling the EVD volume, use n-order B-Spline interpolation.
        int evdScaleSplineOrder = 3;

        // If true, the calculated EVD wont be saved
        bool dontSaveEVD = false;

        /**
         * Load settings from a file.
         *
         * \param fn the filename
         *
         * \return settings
         *
         */
        static Settings fromFile( const std::string& fn );
    };

    namespace __json
    {
        void get( bool& in, std::string name, picojson::value::object::const_iterator v )
        {
            if( ( v->first == name ) && v->second.is< bool >() )
            {
                in = v->second.get< bool >();
                LogD << "Setting - bool - " << name << " : " << std::to_string( in ) << LogEnd;
            }
            else if( v->first == name )
            {
                throw std::runtime_error( std::string( "Parsing error in settings file. Error: " ) + name +
                                          " needs to be a boolean." );
            }
        }

        void get( Real& in, std::string name, picojson::value::object::const_iterator v )
        {
            if( ( v->first == name ) && v->second.is< double >() )
            {
                in = v->second.get< double >();
                LogD << "Setting - double - " << name << " : " << std::to_string( in ) << LogEnd;
            }
            else if( v->first == name )
            {
                throw std::runtime_error( std::string( "Parsing error in settings file. Error: " ) + name +
                                          " needs to be a number." );
            }
        }

        void get( int& in, std::string name, picojson::value::object::const_iterator v )
        {
            if( ( v->first == name ) && v->second.is< double >() )
            {
                in = v->second.get< double >();
                LogD << "Setting - int - " << name << " : " << std::to_string( in ) << LogEnd;
            }
            else if( v->first == name )
            {
                throw std::runtime_error( std::string( "Parsing error in settings file. Error: " ) + name +
                                          " needs to be a number." );
            }
        }
    } // namespace __json

    Settings Settings::fromFile( const std::string& fn )
    {
        std::ifstream t( fn );
        if( !t.is_open() )
        {
            throw std::runtime_error( std::string( "Cannot read settings file." ) );
        }
        std::string str( ( std::istreambuf_iterator< char >( t ) ), std::istreambuf_iterator< char >() );

        picojson::value v;
        auto err = picojson::parse( v, str );
        if( !err.empty() )
        {
            throw std::runtime_error( std::string( "Parsing error in settings file. Error: " ) + err );
        }

        Settings result;

        if( !v.is< picojson::object >() )
        {
            throw std::runtime_error( std::string( "Parsing error in settings file. Error: root is not an object." ) );
        }

        const auto& obj = v.get< picojson::object >();
        for( auto i = obj.begin(); i != obj.end(); ++i )
        {
            __json::get( result.volumeSize, "volumeSize", i );
            __json::get( result.pixelSize, "pixelSize", i );
            __json::get( result.networkIsVoxelSpace, "networkIsVoxelSpace", i );
            __json::get( result.cylindric, "cylindric", i );
            __json::get( result.preCroppedNetwork, "preCroppedNetwork", i );
            __json::get( result.networkCutXY, "networkCutXY", i );
            __json::get( result.networkCutZ, "networkCutZ", i );
            __json::get( result.volumeCutXYVVF, "volumeCutXYVVF", i );
            __json::get( result.volumeCutZVVF, "volumeCutZVVF", i );
            __json::get( result.volumeCutXYEVD, "volumeCutXYEVD", i );
            __json::get( result.volumeCutZEVD, "volumeCutZEVD", i );
            __json::get( result.restrictDegree, "restrictDegree", i );
            __json::get( result.mapLargerToMax, "mapLargerToMax", i );
            __json::get( result.maxBranchpointDegree, "maxBranchpointDegree", i );
            __json::get( result.mergeCloseBranchPoints, "mergeCloseBranchPoints", i );
            __json::get( result.maxMergeDistance, "maxMergeDistance", i );
            __json::get( result.diameterPerLengthMin, "diameterPerLengthMin", i );
            __json::get( result.diameterPerLengthMax, "diameterPerLengthMax", i );
            __json::get( result.diameterPerLengthNumBins, "diameterPerLengthNumBins", i );
            __json::get( result.minEVD, "minEVD", i );
            __json::get( result.useMinEVD, "useMinEVD", i );
            __json::get( result.maxEVD, "maxEVD", i );
            __json::get( result.useMaxEVD, "useMaxEVD", i );
            __json::get( result.evdHistMin, "evdHistMin", i );
            __json::get( result.evdHistMax, "evdHistMax", i );
            __json::get( result.evdNumBins, "evdNumBins", i );
            __json::get( result.evdVolumeSize, "evdVolumeSize", i );
            __json::get( result.evdScaleSplineOrder, "evdScaleSplineOrder", i );
            __json::get( result.evdScaleThreshold, "evdScaleThreshold", i );
            __json::get( result.dontSaveEVD, "dontSaveEVD", i );
        }

        if( !result.cylindric )
        {
            LogI << "Force-Disabled cutting of data. Not allowed for non-cylindric data." << LogEnd;
            result.networkCutZ = 100.0;
            result.networkCutXY = 100.0;
            result.volumeCutXYVVF = 100.0;
            result.volumeCutZVVF = 100.0;
            result.volumeCutXYEVD = 100.0;
            result.volumeCutZEVD = 100.0;
        }

        return result;
    }
}
