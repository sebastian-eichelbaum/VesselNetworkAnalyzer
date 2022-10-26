//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_IO_H
#define NOGO_IO_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Data.h"
#include "Math.h"
#include "Range.h"
#include "Types.h"
#include "Util.h"

#include "Logger.h"
#define LogTag "nogo/IO"

namespace nogo
{
    /**
     * Write the vessel data to a set of python arrays.
     *
     * \param vessels the vessel data
     * \param filename the filename to write to.
     */
    void writePyArrays(const std::string filename, const nogo::VesselNetwork& vessels)
    {
        // Debug Output to file:
        std::ofstream out(filename);

        if (!out.is_open())
        {
            LogE << "Could not open " << filename << " for writing." << LogEnd;
            return;
        }

        // Start printing points
        out << std::setprecision(15);

        out << "# Points. Index referenced by lines. Each point is an array of [x,y,z,radius,uses]. Uses is the amount "
               "of connected edges."
            << std::endl;
        out << "vesselNetworkPoints = [" << std::endl;
        for (auto pointIdx : range(vessels.points.size()))
        {
            out << "  [" << vessels.points[pointIdx][0] << "," << vessels.points[pointIdx][1] << ","
                << vessels.points[pointIdx][2] << "," << vessels.pointRadii[pointIdx] << ","
                << +vessels.pointUsages[pointIdx] << "]," << std::endl;
        }
        out << "]" << std::endl << std::endl;

        out << "# Lines. A line is a set of indices into the point data." << std::endl;
        out << "vesselNetworkLines = [" << std::endl;
        for (const auto& line : vessels.connectedLines)
        {
            out << "  [";
            for (const auto& pIdx : line)
            {
                out << pIdx << ",";
            }
            out << "]," << std::endl;
        }
        out << "]" << std::endl;

        out.close();
    }

    void write(const std::string baseDir, const nogo::VesselNetwork& vessels)
    {
        // forward
        auto filename = baseDir + "VesselNetwork.py";
        LogI << "Writing vessel segments to \"" << filename << "\"." << LogEnd;
        writePyArrays(filename, vessels);
    }

    /**
     * Write a set of points and index-based line-strips to a STLD file.
     *
     * \tparam PointContainerType the associative container used for points
     * \tparam LineStripContainer the associative container storing the line strips as array.
     *
     * \param filename where to write.
     * \param points point list
     * \param lineStrips list of line-strip index lists
     * \param segmentData segment information
     * \param branchPointData branch point info
     * \param accept if returning true, the line strip is accepted.
     */
    template < typename PointContainerType, typename LineStripContainer >
    void writeLineStrips(
        const std::string filename, const PointContainerType& points, const LineStripContainer& lineStrips,
        const SegmentData& segmentData, const BranchPointData& branchPointData,
        std::function< bool(const VesselNetwork::IndexType&) > accept =
            [](const VesselNetwork::IndexType& /* index */) { return true; })
    {
        // Debug Output to file:
        std::ofstream out(filename);

        out << "# This is a/the text-based lines format." << std::endl;
        out << "# " << std::endl;
        out << "# Format: TOKEN VALUE_0 VALUE_1 ... VALUE_n" << std::endl;
        out << "# * Valid tokens are:" << std::endl;
        out << "# -> P X Y Z   -- this is a point with X Y Z coordinates as space separated list." << std::endl;
        out << "# -> L Index_0 Index_1 Index_... -- this is a line-strip consisting of the denoted point indices."
            << std::endl;
        out << "# " << std::endl;
        out << "# * Besides these mandatory tokens, there can be:" << std::endl;
        out << "# -> PA Num Value -- Point attribute. Num defines the number of the attribute (not the point index), "
               "starts at 0."
            << std::endl;
        out << "# -> LA Num Value -- Line attribute. Num defines the number of the attribute (not the line index), "
               "starts at 0."
            << std::endl;
        out << "# -> Name Token String -- Name of a token. Tokens as above, L, P, PA Num, LA Num, ..." << std::endl;
        out << "# -> Description Token String -- Description of a token. Tokens as above, L, P, PA Num, LA Num, ..."
            << std::endl;
        out << "# Now, here is the data you are waiting for ..." << std::endl;
        out << "# " << std::endl;
        out << "# " << std::endl;
        out << "# Names: " << std::endl;
        out << "Name PA 0 "
            << "isBranchPoint" << std::endl;
        out << "Description PA 0 "
            << "An attribute describing whether the point is a BranchPoint ( = 1 ) or not." << std::endl;
        out << "Name PA 1"
            << " isCapillary" << std::endl;
        out << "Description PA 1 "
            << "An attribute describing whether the point is a capillary ( = 1 ) or not." << std::endl;
        out << "Name PA 2"
            << " Degree" << std::endl;
        out << "Description PA 2 "
            << "The branchpoint degree. Non-branchpoints have this set to 0." << std::endl;

        out << "Name LA 0 "
            << "Radius" << std::endl;
        out << "Description LA 0 "
            << "The radius of the segment." << std::endl;
        out << "Name LA 1 "
            << "isCapillary" << std::endl;
        out << "Description LA 1 "
            << "An attribute describing whether the line is a capillary ( = 1 ) or not." << std::endl;
        out << "# " << std::endl;
        out << "# " << std::endl;
        out << "# Data: " << std::endl;

        // Start printing points
        out << std::setprecision(100);
        for (auto pointIdx : range(points.size()))
        {
            out << "P " << points[pointIdx][0] << " " << points[pointIdx][1] << " " << points[pointIdx][2] << std::endl;

            if (branchPointData.idxMap.count(pointIdx))
            {
                auto realIdx = branchPointData.idxMap.at(pointIdx);
                out << "PA 0 " << 1 << std::endl;
                out << "PA 1 " << branchPointData.isCapillary[realIdx] << std::endl;
                out << "PA 2 " << branchPointData.degrees[realIdx] << std::endl;
            }
            else
            {
                out << "PA 0 " << 0 << std::endl;
                out << "PA 1 " << 0 << std::endl;
                out << "PA 2 " << 0 << std::endl;
            }
        }

        // Undo modifications for floats
        out << std::resetiosflags(std::ios_base::basefield | std::ios_base::floatfield | std::ios_base::adjustfield);

        // print all strips now
        for (auto stripIdx : range(lineStrips.size()))
        {
            if (!accept(stripIdx))
            {
                continue;
            }

            out << "L ";
            for (auto pointIdx : lineStrips[stripIdx])
            {
                out << static_cast< size_t >(pointIdx) << " "; // NOTE: the static cast ensures that we do not print
                                                               // chars if the index type might some uint8 or something
                                                               // (which is unlikely but might happen though).
            }
            out << std::endl;
            if (segmentData.radiiSingle.count(stripIdx) != 0)
            {
                out << "LA 0 " << segmentData.radiiSingle.at(stripIdx) << std::endl;
            }
            else
            {
                out << "LA 0 " << -1 << std::endl;
            }

            if (segmentData.idxMap.count(stripIdx) != 0)
            {
                out << "LA 1 " << segmentData.isCapillary[segmentData.idxMap.at(stripIdx)] << std::endl;
            }
            else
            {
                out << "LA 1 " << 0 << std::endl;
            }
        }
    }

    /**
     * Write the given container to a file, using the provided keys.
     *
     * \tparam ContainerType the type of container to print. Basically this can be anything that is iteratable and whose
     * contents provide a operator. \tparam KeyContainerType the type of container to print (as keys). Basically this
     * can be anything that is iteratable and whose contents provide a << operator. \param filename the file where to
     * write to. \param container the container to write \param keyContainer the key-container to write. Needs to have
     * the same size as container, or is empty (no keys written) \param prefix a prefix per value. Empty by default.
     * \param separator a separator per value. Typically this is the value delimiter.
     * \param accept a functor returning false if the element should NOT be written. Useful for selective write.
     */

    template < typename KeyContainerType, typename ContainerType >
    void writeKeyValueList(
        const std::string filename, const KeyContainerType& keyContainer, const ContainerType& container,
        std::function< bool(const VesselNetwork::IndexType&) > accept =
            [](const VesselNetwork::IndexType& /* index */) { return true; },
        const std::string separator = "\n", const std::string prefix = "")
    {
        if ((keyContainer.size() != 0) && (keyContainer.size() != container.size()))
        {
            throw std::runtime_error("Key container needs to have the same size as value container, or no entry.");
        }

        // Open Ascii file
        std::ofstream out(filename);

        // Ensure proper precision
        out << std::setprecision(100);

        // This is completely independent of the used container type
        bool first = true;
        size_t count = 0;
        size_t ignored = 0;
        bool useKeys = (keyContainer.size() != 0);
        for (auto idx : range(container.size()))
        {
            if (accept(idx))
            {
                // if this is the first element of the stream, do not print the separator
                if (!first)
                {
                    out << separator;
                }
                else
                {
                    first = false;
                }

                if (useKeys)
                {
                    out << prefix << keyContainer[idx] << "   " << container[idx];
                }
                else
                {
                    out << prefix << container[idx];
                }
                count++;
            }
            else
            {
                ignored++;
            }
        }
        LogD << "Written " << count << " values to \"" << filename << "\". Ignored " << ignored << " values." << LogEnd;
    }

    /**
     * Write the given container to a file.
     *
     * \tparam ContainerType the type of container to print. Basically this can be anything that is iteratable and whose
     * contents provide a << operator. \param filename the file where to write to. \param container the container to
     * write \param prefix a prefix per value. Empty by default. \param separator a separator per value. Typically this
     * is the value delimiter. \param accept a functor returning false if the element should NOT be written. Useful for
     * selective write.
     */
    template < typename ContainerType >
    void writeValueList(
        const std::string filename, const ContainerType& container,
        std::function< bool(const VesselNetwork::IndexType&) > accept =
            [](const VesselNetwork::IndexType& /* index */) { return true; },
        const std::string separator = "\n", const std::string prefix = "")
    {
        writeKeyValueList(filename, std::vector< int >(), container, accept, separator, prefix);
    }

    /**
     * Write segment data to several files. File name represents the specific data "meaning".
     *
     * \param data the data to write
     * \param baseDir where to write the files.
     */
    void write(const std::string baseDir, const SegmentData& data)
    {
        auto filenameBase = baseDir + "Segments.";
        LogI << "Writing vessel segment information to \"" << filenameBase << "..\"." << LogEnd;

        // for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "lengths.capillary.csv", data.lengths,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "lengths.non-capillary.csv", data.lengths,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "lengths.csv", data.lengths);

        // for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "tortuosity.capillary.csv", data.tortuosity,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "tortuosity.non-capillary.csv", data.tortuosity,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "tortuosity.csv", data.tortuosity);

        // for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "radius.capillary.csv", data.radius,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "radius.non-capillary.csv", data.radius,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "radius.csv", data.radius);

        writeValueList(filenameBase + "angleAdjacency.capillary.csv", data.anglesAdjacency.caps);
        writeValueList(filenameBase + "angleAdjacency.non-capillary.csv", data.anglesAdjacency.noncaps);
        writeValueList(filenameBase + "angleAdjacency.csv", data.anglesAdjacency.all);

        // for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "diameter.capillary.csv", data.diameter,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "diameter.non-capillary.csv", data.diameter,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "diameter.csv", data.diameter);

        // for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "volume.capillary.csv", data.volume,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "volume.non-capillary.csv", data.volume,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "volume.csv", data.volume);

        // density
        writeValueList(filenameBase + "density.capillary.csv", std::vector< Real >(1, data.segmentDensityCapillary));
        writeValueList(filenameBase + "density.non-capillary.csv",
                       std::vector< Real >(1, data.segmentDensityNonCapillary));
        writeValueList(filenameBase + "density.csv", std::vector< Real >(1, data.segmentDensity));

        // num segments
        writeValueList(filenameBase + "count.capillary.csv", std::vector< Real >(1, data.numCapillarySegments));
        writeValueList(filenameBase + "count.non-capillary.csv", std::vector< Real >(1, data.numNonCapillarySegments));
        writeValueList(filenameBase + "count.csv", std::vector< Real >(1, data.numSegments));

        // Bins
        writeKeyValueList(filenameBase + "diameter_per_length_hist.csv", data.meanDiameterPerLengthBinCenters,
                          data.meanDiameterPerLength);
    }

    /**
     * Write segment data to several files. File name represents the specific data "meaning".
     *
     * \param data the data to write
     * \param baseDir where to write the files.
     */
    void write(const std::string baseDir, const BranchPointData& data)
    {
        auto filenameBase = baseDir + "BranchPoints.";
        LogI << "Writing vessel branch point information to \"" << filenameBase << "..\"." << LogEnd;

        // Degrees - for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "degrees.capillary.csv", data.degrees,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "degrees.non-capillary.csv", data.degrees,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "degrees.csv", data.degrees);

        // Degrees as means - for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "degrees_mean.capillary.csv", std::vector< Real >(1, data.degreeMeanCapillary));
        writeValueList(filenameBase + "degrees_mean.non-capillary.csv",
                       std::vector< Real >(1, data.degreeMeanNonCapillary));
        writeValueList(filenameBase + "degrees_mean.csv", std::vector< Real >(1, data.degreeMean));

        // Diameter - for all three cases: capillary, non-capillary, both
        writeValueList(filenameBase + "diameter.capillary.csv", data.diameter,
                       [&](auto index) { return data.isCapillary[index]; });
        writeValueList(filenameBase + "diameter.non-capillary.csv", data.diameter,
                       [&](auto index) { return !data.isCapillary[index]; });
        writeValueList(filenameBase + "diameter.csv", data.diameter);

        // Densities
        writeValueList(filenameBase + "density.capillary.csv",
                       std::vector< Real >(1, data.branchPointDensityCapillary));
        writeValueList(filenameBase + "density.non-capillary.csv",
                       std::vector< Real >(1, data.branchPointDensityNonCapillary));
        writeValueList(filenameBase + "density.csv", std::vector< Real >(1, data.branchPointDensity));

        writeValueList(filenameBase + "density.degree3.capillary.csv",
                       std::vector< Real >(1, data.branchPointDegree3DensityCapillary));
        writeValueList(filenameBase + "density.degree3.non-capillary.csv",
                       std::vector< Real >(1, data.branchPointDegree3DensityNonCapillary));
        writeValueList(filenameBase + "density.degree3.csv", std::vector< Real >(1, data.branchPointDegree3Density));

        writeValueList(filenameBase + "density.degree4.capillary.csv",
                       std::vector< Real >(1, data.branchPointDegree4DensityCapillary));
        writeValueList(filenameBase + "density.degree4.non-capillary.csv",
                       std::vector< Real >(1, data.branchPointDegree4DensityNonCapillary));
        writeValueList(filenameBase + "density.degree4.csv", std::vector< Real >(1, data.branchPointDegree4Density));
    }

    /**
     * Write volume data to several files. File name represents the specific data "meaning".
     *
     * \param data the data to write
     * \param baseDir where to write the files.
     */
    void write(const std::string baseDir, const VolumeData& data)
    {
        auto filenameBase = baseDir + "Volume.";
        LogI << "Writing vessel volume information to \"" << filenameBase << "..\"." << LogEnd;

        // Densities
        writeValueList(filenameBase + "vvf_local.capillary.csv",
                       std::vector< Real >(1, data.vesselVolumeFractionLocalCapillary));
        writeValueList(filenameBase + "vvf_local.non-capillary.csv",
                       std::vector< Real >(1, data.vesselVolumeFractionLocalNonCapillary));
        writeValueList(filenameBase + "vvf_local.csv", std::vector< Real >(1, data.vesselVolumeFractionLocal));

        writeValueList(filenameBase + "vvf.csv", std::vector< Real >(1, data.vesselVolumeFraction));
    }

    /**
     * Write volume evd data to several files. File name represents the specific data "meaning".
     *
     * \param data the data to write
     * \param baseDir where to write the files.
     */
    void write(const std::string baseDir, const EVDData& data)
    {
        auto filenameBase = baseDir + "Volume.";
        LogI << "Writing EVD volume information to \"" << filenameBase << "..\"." << LogEnd;

        // Histogram
        writeKeyValueList(filenameBase + "evd_hist.csv", data.histBinsCenters, data.discretized);
        writeKeyValueList(filenameBase + "evd_relative_hist.csv", data.histBinsCenters, data.discretizedRelative);
        writeValueList(filenameBase + "evd_mean.csv", std::vector< Real >(1, data.mean));
    }

    /**
     * Write segment data to several files. File name represents the specific data "meaning".
     *
     * \param segments the segments of the vessel network
     * \param points the raw network point data loaded.
     * \param baseDir where to write the files.
     * \param segmentData segment information
     * \param branchPointData branch point info
     * \param filenamePrefix a prefix to append to the filename. Can be omitted.
     */
    void write(const std::string baseDir, const VesselSegments& segments, const decltype(VesselNetwork::points)& points,
               const SegmentData& segmentData, const BranchPointData& branchPointData, std::string filenamePrefix = "")
    {
        // forward
        auto filename = baseDir + filenamePrefix + "VesselSegments.stld";
        LogI << "Writing vessel segments to \"" << filename << "\"." << LogEnd;
        writeLineStrips(filename, points, segments.segments, segmentData, branchPointData,
                        [&](auto&& segmentIdx) { return !segments.removed[segmentIdx]; });
    }

    /**
     * Load the network from a VTK PolyData file.
     *
     * \param filename the filename
     * \param scale Scale the network using this factor
     *
     * \throw std::runtime_error if the file is somehow invalid or the file could not be opened.
     *
     * \return the vessel network as represented by the file
     */
    inline SPtr< VesselNetwork > loadVesselNetworkVTK(const std::string& filename, double scale = 1.0)
    {
        // The VTK types to use
        using VTKIndexType = int;
        using VTKValueType = float;
        using VTKPointType = std::array< VTKValueType, 3 >;

        // Load into those structures first, as the target value type in VesselNetwork might be different.
        std::vector< VTKPointType > loadedPoints;
        std::vector< VTKValueType > loadedRadii;
        std::vector< VTKValueType > loadedPointRadii;
        std::vector< VTKIndexType > loadedLines;

        // Open file
        std::ifstream input(filename, std::ios::binary);
        if (!input.good())
        {
            throw std::runtime_error("File \"" + filename + "\" could not be opened for reading.");
        }

        // get length of file:
        input.seekg(0, input.end);
        int length = input.tellg();
        input.seekg(0, input.beg);

        // Read lines until binary blocks are reached
        // NOTE: the first 5 lines of a VTK file are ASCII, regardless of binary or ascii data modes. Read and check
        // validity:
        std::string line;
        std::getline(input, line);
        if (((line != "# vtk DataFile Version 3.0") && (line != "# vtk DataFile Version 4.1")) || !input.good())
        {
            throw std::runtime_error("File \"" + filename + "\", line 1 needs to be vtk version header.");
        }

        // Second line is a name/description -> skip
        std::getline(input, line);
        if (!input.good())
        {
            throw std::runtime_error("File \"" + filename + "\" unexpected EOF.");
        }

        // Third line denotes binary/ascii mode
        std::getline(input, line);
        if (line != "BINARY" || !input.good())
        {
            throw std::runtime_error("File \"" + filename + "\" is not a binary VTK file.");
        }

        // 4th line denotes dataset type -> we require polydata here
        std::getline(input, line);
        if (line != "DATASET POLYDATA" || !input.good())
        {
            throw std::runtime_error("File \"" + filename + "\" is not a PolyData VTK file.");
        }

        LogD << "Header loaded. Binary PolyData file. Continuing to load blocks." << LogEnd;

        // Iterate the remaining file and evaluate the blocks to follow.
        VesselNetwork::IndexType numLines = 0;
        bool hasLines = false;
        bool hasPoints = false;
        bool hasRadii = false;
        bool hasPointRadii = false;
        bool hasPolylines = false;

        while (input.good() && !input.eof())
        {
            LogD << "Position: " << input.tellg() << LogEnd;

            // Read Block Header
            std::getline(input, line);
            auto elems = split(line);
            if (elems.size() < 1)
            {
                throw std::runtime_error("File \"" + filename + "\", at byte " + std::to_string(input.tellg()) +
                                         " contains invalid block descriptor \"" + line + "\".");
            }

            // Is POINTS?
            if (elems[0] == "POINTS")
            {
                if (elems.size() != 3)
                {
                    throw std::runtime_error("File \"" + filename + "\" contains invalid Points block header \"" +
                                             line + "\".");
                }

                LogD << "Points Block found: " << elems[1] << " points of type " << elems[2] << LogEnd;
                if (elems[2] != "float")
                {
                    throw std::runtime_error("File \"" + filename + "\" contains non-float points block.");
                }
                auto numElements = std::atol(elems[1].c_str());

                // Alloc some mem and copy to vector.
                // NOTE: use resize as reserve does not update the internal size counter.
                loadedPoints.resize(numElements);
                input.read(reinterpret_cast< char* >(loadedPoints.data()), numElements * sizeof(VTKPointType));

                hasPoints = true;

                // According to the VTK doc, a block ends with a line break. Skip the break and continue reading the
                // next block.
                input.ignore();
            }
            else if (elems[0] == "LINES") // Is LINES?
            {
                if (elems.size() != 3)
                {
                    throw std::runtime_error("File \"" + filename + "\" contains invalid Lines block header \"" + line +
                                             "\".");
                }

                LogD << "Lines Block found: " << elems[1] << " lines, " << elems[2] << " cell-blocksize." << LogEnd;
                numLines = std::atol(elems[1].c_str());
                auto blockSize = std::atol(elems[2].c_str());

                // Alloc some mem and copy to vector
                // NOTE: use resize as reserve does not update the internal size counter.
                loadedLines.resize(blockSize);
                input.read(reinterpret_cast< char* >(loadedLines.data()), blockSize * sizeof(VTKIndexType));

                hasLines = true;
                // If a line is NOT only index and 2 point indices (=3), assume the data is made up of polylines
                hasPolylines = ((numLines * 3) != blockSize);

                // According to the VTK doc, a block ends with a line break. Skip the break and continue reading the
                // next block.
                input.ignore();
            }
            else if (elems[0] == "CELL_DATA") // Is additional data per cell??
            {
                if (elems.size() != 2)
                {
                    throw std::runtime_error("File \"" + filename + "\" contains invalid Cell Data block header \"" +
                                             line + "\".");
                }

                LogD << "Cell Data Block found: " << elems[1] << " values." << LogEnd;

                // The Cell_Data block contains another header:
                std::getline(input, line);
                if (line != "FIELD FieldData 1" || !input.good())
                {
                    throw std::runtime_error("File \"" + filename +
                                             "\" has a cell data block with unsupported field type: \"" + line + "\"");
                }

                // And the Field itself: Name Index NumValues Type
                std::getline(input, line);
                auto fieldElems = split(line);
                if (fieldElems.size() != 4)
                {
                    throw std::runtime_error("File \"" + filename + "\" uses invalid field data header: \"" + line +
                                             "\"");
                }

                LogD << "Field Data found: " << fieldElems[0] << " with " << fieldElems[2] << " values of type "
                     << fieldElems[3] << LogEnd;
                if (fieldElems[0] != "radius")
                {
                    throw std::runtime_error("File \"" + filename + "\" contains field data that is not \"radius\".");
                }

                if (fieldElems[3] != "float")
                {
                    throw std::runtime_error("File \"" + filename + "\" contains non-float field data block.");
                }
                auto numValues = std::atol(fieldElems[2].c_str());

                // Alloc some mem and copy to vector.
                // NOTE: use resize as reserve does not update the internal size counter.
                loadedRadii.resize(numValues);
                input.read(reinterpret_cast< char* >(loadedRadii.data()), numValues * sizeof(VTKValueType));

                hasRadii = true;

                // According to the VTK doc, a block ends with a line break. Skip the break and continue reading the
                // next block.
                input.ignore();
            }
            else if (elems[0] == "POINT_DATA") // Ignore it parse point data?
            {
                if (elems.size() != 2)
                {
                    throw std::runtime_error("File \"" + filename + "\" contains invalid Point Data block header \"" +
                                             line + "\".");
                }

                LogD << "Point Data Block found: " << elems[1] << " values." << LogEnd;

                // The Point block contains another header:
                std::getline(input, line);
                if (line != "FIELD FieldData 1" || !input.good())
                {
                    throw std::runtime_error("File \"" + filename +
                                             "\" has a point data block with unsupported field type: \"" + line + "\"");
                }

                // And the Field itself: Name Index NumValues Type
                std::getline(input, line);
                auto fieldElems = split(line);
                if (fieldElems.size() != 4)
                {
                    throw std::runtime_error("File \"" + filename + "\" uses invalid field data header: \"" + line +
                                             "\"");
                }

                LogD << "Field Data found: " << fieldElems[0] << " with " << fieldElems[2] << " values of type "
                     << fieldElems[3] << LogEnd;
                if (fieldElems[0] != "radius")
                {
                    throw std::runtime_error("File \"" + filename + "\" contains field data that is not \"radius\".");
                }

                if (fieldElems[3] != "float")
                {
                    throw std::runtime_error("File \"" + filename + "\" contains non-float field data block.");
                }
                auto numValues = std::atol(fieldElems[2].c_str());

                // Alloc some mem and copy to vector.
                // NOTE: use resize as reserve does not update the internal size counter.
                loadedPointRadii.resize(numValues);
                input.read(reinterpret_cast< char* >(loadedPointRadii.data()), numValues * sizeof(VTKValueType));

                hasPointRadii = true;

                // According to the VTK doc, a block ends with a line break. Skip the break and continue reading the
                // next block.
                input.ignore();
            }
            else
            {
                throw std::runtime_error("File \"" + filename + "\" contains unsupported block: \"" + line + "\".");
            }

            // Be nice and also allow one additional char at the end, since some files contain another byte at the end.
            if (input.tellg() >= length - 2)
            {
                break;
            }
        }

        if (!(hasRadii && hasPoints && hasLines))
        {
            throw std::runtime_error("File \"" + filename + "\" does not contain one of the blocks required. POINTS: " +
                                     std::to_string(hasPoints) + " LINES: " + std::to_string(hasLines) +
                                     " RADII: " + std::to_string(hasRadii) + ".");
        }

        if (!(hasPolylines && hasPointRadii))
        {
            throw std::runtime_error("File \"" + filename + "\" contains polylines but no radius per point.");
        }

        // Copy to target structure. Cast types and apply byteSwap.
        auto result = std::make_shared< VesselNetwork >();

        // First, change byte order of the points data, since VTK saves binary files in big-endian order.
        result->points.reserve(loadedPoints.size());
        Vec3 bbMin = {{std::numeric_limits< VesselNetwork::ValueType >::max(),
                       std::numeric_limits< VesselNetwork::ValueType >::max(),
                       std::numeric_limits< VesselNetwork::ValueType >::max()}};
        Vec3 bbMax = {{std::numeric_limits< VesselNetwork::ValueType >::min(),
                       std::numeric_limits< VesselNetwork::ValueType >::min(),
                       std::numeric_limits< VesselNetwork::ValueType >::min()}};
        for (auto point : loadedPoints)
        {
            // *pointIter is a point, consisting out of 3 vtk values -> swap the values of the vector
            auto swappedPoint = byteSwap(point, sizeof(VTKValueType));

            // use this chance to calc the BB
            auto pX = scale * static_cast< VesselNetwork::ValueType >(swappedPoint[0]);
            auto pY = scale * static_cast< VesselNetwork::ValueType >(swappedPoint[1]);
            auto pZ = scale * static_cast< VesselNetwork::ValueType >(swappedPoint[2]);
            bbMax[0] = std::max(pX, bbMax[0]);
            bbMax[1] = std::max(pY, bbMax[1]);
            bbMax[2] = std::max(pZ, bbMax[2]);
            bbMin[0] = std::min(pX, bbMin[0]);
            bbMin[1] = std::min(pY, bbMin[1]);
            bbMin[2] = std::min(pZ, bbMin[2]);

            result->points.push_back({{pX, pY, pZ}});
        }

        result->boundingBox.first = bbMin;
        result->boundingBox.second = bbMax;

        // Next, change byte order of the radius data, since VTK saves binary files in big-endian order.
        if (!hasPolylines)
        {
            result->radii.reserve(loadedRadii.size());
            for (auto radius : loadedRadii)
            {
                result->radii.push_back(scale * static_cast< VesselNetwork::ValueType >(byteSwap(radius)));
            }
        }
        else
        {
            result->radii.reserve(loadedPointRadii.size());
        }

        // We now have the line-strips as a block of run-length encoded indices. For a more comfortable use, we
        // re-arrange it a bit. -> inflate the run-length encoding.
        size_t idx = 0;
        // size_t lIdx = 0;
        result->pointRadii = std::vector< Real >(loadedPoints.size(), Real(0));
        result->pointUsages = decltype(result->pointUsages)(loadedPoints.size(), 0);
        while (loadedLines.size() > idx)
        {
            auto len = byteSwap(loadedLines[idx++]);
            // const auto lineRadius = static_cast< VesselNetwork::ValueType >( byteSwap( loadedRadii[ lIdx++ ] ) );

            result->connectedLines.emplace_back();
            result->connectedLines.back().reserve(len);
            for (const auto pI : range(len - 1))
            {
                auto vI1 = byteSwap(loadedLines[idx++]);
                auto vI2 = byteSwap(loadedLines[idx]);

                // This ensures we have valid indices
                if (static_cast< size_t >(vI1) >= result->points.size())
                {
                    throw std::runtime_error(
                        "File \"" + filename +
                        "\" contains lines with point indices out of range: " + std::to_string(vI1) +
                        " >= " + std::to_string(result->points.size()) + " for " + std::to_string(pI) + ".");
                }
                if (static_cast< size_t >(vI2) >= result->points.size())
                {
                    throw std::runtime_error(
                        "File \"" + filename +
                        "\" contains lines with point indices out of range: " + std::to_string(vI2) +
                        " >= " + std::to_string(result->points.size()) + " for " + std::to_string(pI) + ".");
                }

                result->pointRadii[vI1] =
                    scale * static_cast< VesselNetwork::ValueType >(byteSwap(loadedPointRadii[vI1]));
                result->pointRadii[vI2] =
                    scale * static_cast< VesselNetwork::ValueType >(byteSwap(loadedPointRadii[vI2]));
                result->pointUsages[vI1]++;
                result->pointUsages[vI2]++;

                const auto r1 = result->pointRadii[vI1];
                const auto r2 = result->pointRadii[vI2];
                result->radii.push_back((r1 + r2) / 2.0);
                // result->radii.push_back( r2 );
                // result->radii.push_back( scale * lineRadius );

                result->lines.push_back(std::make_pair(vI1, vI2));
                result->connectedLines.back().push_back(vI1);
                if (pI == (len - 2))
                {
                    // Fore add the last point
                    result->connectedLines.back().push_back(vI2);
                }
            }
            idx++;
        }

        /*
        // We now have the line-strips as a block of run-length encoded indices. For a more comfortable use, we
        // re-arrange it a bit. -> inflate the run-length encoding.
        size_t idx = 0;
        while (loadedLines.size() > idx)
        {
            auto len = byteSwap(loadedLines[idx]);
            if (len != 2)
            {
                throw std::runtime_error("File \"" + filename +
                                         "\" contains line-strips with more than one segment. This is not supported.");
            }

            // Add to lines list
            auto i1 = byteSwap(loadedLines[idx + 1]);
            auto i2 = byteSwap(loadedLines[idx + 2]);
            result->lines.push_back(std::make_pair(i1, i2));

            // This ensures we have valid indices
            if (static_cast< size_t >(i1) >= result->points.size())
            {
                throw std::runtime_error("File \"" + filename + "\" contains lines with point indices out of range: " +
                                         std::to_string(i1) + " >= " + std::to_string(result->points.size()) + ".");
            }
            if (static_cast< size_t >(i2) >= result->points.size())
            {
                throw std::runtime_error("File \"" + filename + "\" contains lines with point indices out of range: " +
                                         std::to_string(i2) + " >= " + std::to_string(result->points.size()) + ".");
            }

            // Next index is line length and the additional slot for the length value itself:
            idx += len + 1;
        }
        */

        LogI << "Loaded VTK PolyData file in binary format. " << length << " bytes loaded." << LogEnd;
        return result;
    }

    /**
     * Load volume file into a volume structure.
     *
     * \param filename file to load
     * \param volume the volume to fill
     *
     * \return the volume
     */
    template < typename ValueType >
    void loadVolume(const std::string& filename, SPtr< Volume< ValueType > > volume)
    {
        // Open the file
        std::ifstream f(filename,
                        std::ios::in | std::ios::ate | std::ios::binary); // NOTE: ios::ate sets pointer to end of file
        if (!f.good())
        {
            throw std::runtime_error("Could not open file \"" + filename + "\" for reading.");
        }

        // Determine the file length
        std::size_t sizeBytes = f.tellg();
        f.seekg(0);

        // How many items?
        std::size_t numItems = sizeBytes / sizeof(ValueType);

        // Load the data
        volume->data.resize(numItems);
        f.read(reinterpret_cast< char* >(&volume->data[0]), sizeBytes);

        LogD << "Loaded " << sizeBytes << " bytes ( = " << numItems << " values) from file." << LogEnd;
        f.close();
    }

    /**
     * Write the given volume to a file.
     *
     * \param filename the file to save to
     * \param data the data to write
     */
    template < typename ValueType >
    void write(const std::string filename, const Volume< ValueType >& volume)
    {
        // Open the file
        std::ofstream f(filename, std::ios::out | std::ios::binary);
        if (!f.good())
        {
            throw std::runtime_error("Could not open file \"" + filename + "\" for writing.");
        }

        auto size = volume.size();
        f.write(reinterpret_cast< const char* >(&volume.data[0]), size * sizeof(ValueType));

        LogD << "Written " << size * sizeof(ValueType) << " bytes ( " << size << " items ) to file." << LogEnd;
        f.close();
    }
} // namespace nogo

#endif // NOGO_IO_H
