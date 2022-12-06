//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#ifndef NOGO_DATA_H
#define NOGO_DATA_H

#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "Math.h"
#include "Types.h"

#include "Logger.h"
#define LogTag "nogo/Data"

namespace nogo
{
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // INPUT data
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template < typename ValueType >
    struct Volume
    {
        /**
         * Type used for indexing values in the data
         */
        using IndexType = int64_t;

        /**
         * he type of values stored
         */
        using value_type = ValueType;

        /**
         * Create a volume of a given size
         *
         * \param sX size in X
         * \param sY size in Y
         * \param sZ size in Z
         */
        Volume(const IndexType& sX, const IndexType& sY, const IndexType& sZ)
            : sizeX(sX), sizeY(sY), sizeZ(sZ),
              boundingBox(
                  std::make_pair(Vec3{{static_cast< Real >(0), static_cast< Real >(0), static_cast< Real >(0)}},
                                 Vec3{{static_cast< Real >(sX), static_cast< Real >(sY), static_cast< Real >(sZ)}}))
        {
            data.resize(size());
        }

        /**
         * The size of the volume
         */
        IndexType sizeX, sizeY, sizeZ;

        /**
         * Calculate the index in the data array for a given coordimate
         *
         * \param x x coord
         * \param y y coord
         * \param z z coord
         *
         * \return memory index in the data
         */
        IndexType index(const IndexType& x, const IndexType& y, const IndexType& z) const
        {
            return (z * sizeY * sizeX) + (y * sizeX) + x;
        }

        /**
         * Checks if the given coordinate is valid.
         *
         * \param x x coord
         * \param y y coord
         * \param z z coord
         *
         * \return true if valid
         */
        bool valid(const IndexType& x, const IndexType& y, const IndexType& z) const
        {
            auto lowerCriterion = (x >= 0) && (y >= 0) && (z >= 0);
            auto upperCriterion = (x < sizeX) && (y < sizeY) && (z < sizeZ);
            return upperCriterion && lowerCriterion;
        }

        /**
         * Size of the volume.
         *
         * \return the size as number of elements.
         */
        IndexType size() const
        {
            return sizeX * sizeY * sizeZ;
        }

        /**
         * The volume itself.
         */
        std::vector< ValueType > data;

        /**
         * The data bounding box in voxel coordinates.
         */
        BoundingBox boundingBox;
    };

    /**
     * Lazy abbreviations
     */
    using MaskVolume = Volume< uint8_t >;
    using DTVolume = Volume< float >;
    using CTVolume = Volume< uint16_t >;

    /**
     * Basic data capsule for vessel network data. General rule: DO NOT write methods that derive data from data in the
     * class. Only provide methods to access data.
     */
    struct VesselNetwork
    {
        /**
         * Type used for indexing in line-strips
         */
        using IndexType = int;

        /**
         * Type used for storing attributes (radius).
         */
        using ValueType = Real;

        /**
         * The points in space. Vector index is point index.
         */
        std::vector< Vec3 > points;

        /**
         * The radii associated with each line. A line index in m_lines matches its radius index.
         */
        std::vector< ValueType > radii;

        /**
         * The radii associated with each point. A point index matches its radius index.
         */
        std::vector< ValueType > pointRadii;

        /**
         * The lines in the graph as pair of node indices. Multiple lines might make up a edge (between two branch
         * points).
         */
        std::vector< std::pair< IndexType, IndexType > > lines;

        /**
         * Connected lines. A line is a set of point indices.
         */
        std::vector< std::vector< IndexType > > connectedLines;

        /**
         * Counts the use of each point.
         */
        std::vector< uint8_t > pointUsages;

        /**
         * The bounding box of the data.
         */
        BoundingBox boundingBox;

        /**
         * Use an explizit volume for that network. If smaller than zero, the AABB volume will be used instead.
         */
        Real volumeOverride = static_cast< Real >(-1);
    };

    /**
     * Derived from the input data, this describes the actual vessel segments. The input data usually is only a bunch of
     * lines. The real adjacency information is given implicitly only. This struct contains adjacency explicitly.
     */
    struct VesselSegments
    {
        /**
         * The structure for storing index lists. This represents adjacent lines, connected by a shared point.
         */
        using Segment = std::list< VesselNetwork::IndexType >;

        /**
         * Type used to store multiple radii of a segment.
         */
        using SegmentRadii = std::map< Real, VesselNetwork::IndexType >;

        //! The radii associated with a single point.
        using PointRadii = std::vector< VesselNetwork::ValueType >;

        /**
         * A copy of the point list, as the network reconstruction algorithm modifies the data.
         */
        decltype(VesselNetwork::points) points;

        /**
         * Point degrees.
         */
        std::vector< size_t > degrees;

        /**
         * The list of segments, representing ALL segments.
         */
        std::vector< Segment > segments;

        /**
         * Removed segments.
         */
        std::vector< bool > removed;

        /**
         * The radii per segment. It might be the case that a segment has multiple radii. Count the several radii.
         */
        std::vector< SegmentRadii > radii;

        /**
         * The radii associated with each point. A point index matches its radius index. This is not necessarily the
         * same as VesselNetwork::pointRadii - it is the list of the radii of each line adjacent to the given point
         * index.
         */
        std::vector< PointRadii > pointRadii;

        /**
         * Length of a segment in euclidean space.
         */
        std::vector< Real > length;

        /**
         * Capillary volume derived from lines data.
         */
        Real volumeCapillaries;

        /**
         * Non-Capillary volume derived from lines data.
         */
        Real volumeNonCapillaries;
    };

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // OUTPUT data
    //
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * This struct contains the local segment data like lengths and so on. It is just a container for several values
     * associated with vessel segments. It also contains some kind of redundant information. I.e. the radius AND the
     * diameter. The purpose is to ensure that the scripts used to generate the diagrams do not need to calculate
     * anything. They just take these arrays and draw diagrams.
     */
    struct SegmentData
    {
        /**
         * Map between the real segment index and the index in the containers in here.
         */
        std::map< VesselNetwork::IndexType, VesselNetwork::IndexType > idxMap;

        /**
         * The length of the segment. The index in the vector relates to the segment ID.
         */
        std::vector< Real > lengths;

        /**
         * The direct distance between start and endpoint of the segment. The index in the vector relates to the segment
         * ID. Although not directly needed for the requested figures, it might get handy some time.
         */
        std::vector< Real > directLengths;

        /**
         * The relation between direct length to length. The index in the vector relates to the segment ID.
         */
        std::vector< Real > tortuosity;

        /**
         * Segment radius. The index in the vector relates to the segment ID.
         */
        std::vector< Real > radius;

        /**
         * Segment diameter. Its always 2.0 * radius at the given index. The index in the vector relates to the segment
         * ID.
         */
        std::vector< Real > diameter;

        /**
         * Cylindric volume of a segment.
         */
        std::vector< Real > volume;

        struct CapNonCapVector
        {
            std::vector< Real > all;
            std::vector< Real > caps;
            std::vector< Real > noncaps;

            void addByDiameter(Real /*d*/, Real angle, bool isCapillary)
            {
                all.push_back(angle);
                if (isCapillary)
                {
                    caps.push_back(angle);
                }
                else
                {
                    noncaps.push_back(angle);
                }
            }

            void reserve(size_t n)
            {
                all.reserve(n);
                caps.reserve(n);
                noncaps.reserve(n);
            }
        };

        //! Angles of vectors being adjacent to each other at a branch point
        CapNonCapVector anglesAdjacency;

        /**
         * Is this segment a capillary according a defined predicate? The index in the vector relates to the segment ID.
         */
        std::vector< bool > isCapillary;

        /**
         * A merged radius per segment.
         */
        std::map< VesselNetwork::IndexType, Real > radiiSingle;

        /**
         * The mean diameter per length-span.
         */
        std::vector< Real > meanDiameterPerLength;

        /**
         * The centers of the length-spans.
         */
        std::vector< Real > meanDiameterPerLengthBinCenters;

        /**
         * Overall segment density
         */
        Real segmentDensity;

        /**
         * Segment density of capillaries
         */
        Real segmentDensityCapillary;

        /**
         * Non-Cap segment density.
         */
        Real segmentDensityNonCapillary;

        //! Number of segments
        size_t numSegments;
        size_t numCapillarySegments;
        size_t numNonCapillarySegments;
    };

    /**
     * This struct contains the local branch point data like degree, diameters and so on. It is only a container for all
     * the data needed.
     */
    struct BranchPointData
    {
        /**
         * Map between the point index list in VesselSegments::points and the indices used here.
         */
        std::map< VesselNetwork::IndexType, VesselNetwork::IndexType > idxMap;

        /**
         * The degrees of all the bifurcations ( branch-points ). Index in the vector is the branch-point ID.
         */
        std::vector< Real > degrees;

        /**
         * Degree mean
         */
        Real degreeMean;

        /**
         * Degree mean
         */
        Real degreeMeanCapillary;

        /**
         * Degree mean
         */
        Real degreeMeanNonCapillary;

        /**
         * The diameter of each branch point.  Index in the vector is the branch-point ID.
         */
        std::vector< Real > diameter;

        /**
         * Is this segment a capillary according a defined predicate? Index in the vector is the branch-point ID.
         */
        std::vector< bool > isCapillary;

        /**
         * Density of branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDensity;

        /**
         * Density of degree 3 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree3Density;

        /**
         * Density of degree 4 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree4Density;

        /**
         * Density of branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDensityCapillary;

        /**
         * Density of degree 3 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree3DensityCapillary;

        /**
         * Density of degree 4 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree4DensityCapillary;

        /**
         * Density of branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDensityNonCapillary;

        /**
         * Density of degree 3 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree3DensityNonCapillary;

        /**
         * Density of degree 4 branch-points in mm^-3. Relative to the cylinder volume.
         */
        Real branchPointDegree4DensityNonCapillary;
    };

    /**
     * Store derived volume information.
     */
    struct VolumeData
    {
        /**
         * The fraction between vessel and non-vessel voxels.
         */
        Real vesselVolumeFraction;

        /**
         * The fraction between vessel and non-vessel voxels. In the cylinder of the network.
         */
        Real vesselVolumeFractionLocal;

        /**
         * The fraction between vessel and non-vessel voxels. For capillaries only.
         */
        Real vesselVolumeFractionLocalCapillary;

        /**
         * The fraction between vessel and non-vessel voxels. For non-capillaries only.
         */
        Real vesselVolumeFractionLocalNonCapillary;
    };

    /**
     * Store derived EVD data. This actually is a histogram of the distance transformed volume.
     */
    struct EVDData
    {
        /**
         * Max EVD
         */
        Real max;

        /**
         * Mean EVD
         */
        Real mean;

        /**
         * Relative histogram of all EVD data
         */
        std::vector< Real > discretizedRelative;

        /**
         * Absolute histogram of all EVD data.
         */
        std::vector< size_t > discretized;

        /**
         * histogram minimum.
         */
        Real histMin;

        /**
         * Histogram max
         */
        Real histMax;

        /**
         * Num of bins
         */
        size_t histNumBins;

        /**
         * Bin centers
         */
        std::vector< Real > histBinsCenters;
    };
} // namespace nogo

#endif // NOGO_DATA_H
