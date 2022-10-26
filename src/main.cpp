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
#include <string>

#include "Analyze.h"
#include "Data.h"
#include "DistanceTransform.h"
#include "IO.h"
#include "Settings.h"

#include "Logger.h"
#define LogTag "main"

namespace
{
    /**
     * Print welcome message and copyright information.
     */
    void printHello()
    {
        std::cout << "Vessel Network Analyzer" << std::endl
                  << "Copyright (C) 2015-2021 Dr. Sebastian Eichelbaum (http://www.nemtics.com)" << std::endl
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
        std::cout << "Usage:" << std::endl
                  << std::endl
                  << "analyzer OUTPUT_PATH NETWORK_FILE [VOLUME_FILE] [VOLUMEDT_FILE] [SETTINGS_FILE]" << std::endl
                  << std::endl
                  << "  OUTPUT_PATH   - the relative or absolute path where to write results. Needs to exist."
                  << std::endl
                  << "  NETWORK_FILE  - the relative or absolute path to the network file." << std::endl
                  << "  VOLUME_FILE   - optional - the relative or absolute path to the volume file. To skip, specify "
                     "an empty string via \"\"."
                  << std::endl
                  << "  VOLUMEDT_FILE - optional - the relative or absolute path to the distance transformed volume "
                     "file. To skip, specify an empty string via \"\"."
                  << std::endl
                  << "  SETTINGS_FILE - optional - the config file to load. To skip, specify an empty string via \"\"."
                  << std::endl;
    }
} // namespace

/**
 * The program's main function and entry point. Start reading here.
 *
 * \param argc number of arguments given on command line
 * \param argv argument list
 *
 * \return 0 if everything was fine. Return something else to denote an error.
 */
int main(int argc, char** argv)
{
    printHello();

    // Parameter validity
    if ((argc != 6))
    {
        printUsage();
        return -1;
    }

    LogI << "Startup successful. Starting ..." << LogEnd;
    LogD << "Parameter Summary:" << LogEnd;
    for (auto argumentIndex : nogo::range(0, argc))
    {
        LogD << argumentIndex << ": " << argv[argumentIndex] << LogEnd;
    }

    // Catch any kind of exception
    try
    {
        // 1 - Loading the data
        auto networkFile = std::string(argv[2]);
        auto outputDir = std::string(argv[1]);

        // Use EVD calculations if the DT file was specified
        auto volumeDTFile = std::string(argv[4]);
        auto useEVD = !volumeDTFile.empty();

        // Use VVF calculations if the mask file was specified
        auto volumeFile = std::string(argv[3]);
        auto useVVF = !volumeFile.empty();

        // Used settings
        auto settingsFile = std::string(argv[5]);

        if (outputDir.empty())
        {
            throw std::runtime_error("Output path is not specified.");
        }
        if ((outputDir.back() != '/') && (outputDir.back() != '\\'))
        {
            outputDir += "/";
        }

        // Load settings
        nogo::Settings settings;
        if (!settingsFile.empty())
        {
            LogI << "Try loading settings from \"" << settingsFile << "\"." << LogEnd;
            try
            {
                settings = nogo::Settings::fromFile(settingsFile);
            }
            catch (std::exception& e)
            {
                LogE << "Loading settings failed. Reason: " << e.what() << LogEnd;
                return -1;
            }
        }

        // Start off by loading the network for topological analysis
        LogI << "Try loading vessel network from \"" << networkFile << "\"." << LogEnd;
        auto vessels = nogo::loadVesselNetworkVTK(networkFile, settings.networkIsVoxelSpace ? settings.pixelSize : 1.0);

        // Write as python arrays?
        // LogI << "Save vessel network to \"" << networkFile << "\"." << LogEnd;
        // nogo::write( outputDir, *vessels );

        nogo::SPtr< nogo::MaskVolume > volume;
        if (useVVF)
        {
            LogI << "Try loading volume from \"" << volumeFile << "\" as " << settings.volumeSize << "^3." << LogEnd;
            try
            {
                volume =
                    std::make_shared< nogo::MaskVolume >(settings.volumeSize, settings.volumeSize, settings.volumeSize);
                nogo::loadVolume(volumeFile, volume);
            }
            catch (std::exception& e)
            {
                LogE << "Loading mask volume failed. Reason: " << e.what() << LogEnd;
                return -1;
            }
        }

        // Calculate the DT volume?
        nogo::SPtr< nogo::DTVolume > dtVolume;
        bool dtCalced = false;
        if (useEVD)
        {
            try
            {
                LogI << "Try loading dt volume from \"" << volumeDTFile << "\"." << LogEnd;
                auto vsize = (settings.evdVolumeSize > 0) ? std::min(settings.volumeSize, settings.evdVolumeSize)
                                                          : settings.volumeSize;
                dtVolume = std::make_shared< nogo::DTVolume >(vsize, vsize, vsize);
                nogo::loadVolume(volumeDTFile, dtVolume);
            }
            catch (std::exception& e)
            {
                dtVolume.reset();

                // The volume can only be calculated if the mask volume was specified too.
                if (useVVF)
                {
                    LogI << "Distance Transformed volume could not be loaded (Reason: " << e.what()
                         << "). Calculating. This might take some serious time!" << LogEnd;
                    dtVolume = nogo::dt(*volume, settings.evdVolumeSize, settings.evdScaleSplineOrder,
                                        settings.evdScaleThreshold);
                    dtCalced = true;
                }
                else
                {
                    LogE << "Cannot calculate EVD without a mask volume. Specify mask volume file to force EVD "
                            "calculation."
                         << LogEnd;
                    return -1;
                }
            }
            if (dtCalced && !settings.dontSaveEVD)
            {
                LogI << "Calculated distance transformed volume. Writing to \"" << volumeDTFile << "\"." << LogEnd;
                nogo::write(volumeDTFile, *dtVolume);
            }
        }

        // 2 - Process
        nogo::analyze(outputDir, *vessels, volume, dtVolume, settings);

        LogI << "Analysis done. Quit." << LogEnd;
    }
    catch (std::exception& e)
    {
        LogE << "Exception caught: " << e.what() << LogEnd;
        return -1;
    }

    LogI << "Bye." << LogEnd;

    return 0;
}
