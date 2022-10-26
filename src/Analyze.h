//----------------------------------------------------------------------------------------
//
// Project: Analyzer
//
// Copyright (C) 2015 Sebastian Eichelbaum (http://www.nemtics.com)
//
// You should have received a copy of the License along with this program.
//
//----------------------------------------------------------------------------------------

#include <algorithm>
#include <list>
#include <map>
#include <numeric>
#include <vector>

// For debugging to write files:
#include <fstream>

#include "Data.h"
#include "IO.h"
#include "Math.h"
#include "Range.h"
#include "Settings.h"
#include "Types.h"

#include "Logger.h"
#define LogTag "nogo/Analyze"

namespace nogo
{
    using BB = decltype(VesselNetwork::boundingBox);

    constexpr auto pi = 3.14159265359;

    // A predicate denoting whether the radius denotes a capillary or not.
    bool isCapillary(const VesselNetwork::ValueType& radius)
    {
        return (radius < 3.5);
    }

    /**
     * Cut the given bounding box using a given percentage. The cut operation is done along the centerline of the
     * cylinder represented by the bb.
     *
     * \param bb the bb to cut
     * \param cutPerc the percentage
     * \param cutPercZ the percentage to cut in height
     *
     * \return the cut bb
     */
    BB cutBoundingBox(const BoundingBox& bb, Real cutPerc, Real cutPercZ = 100.0)
    {
        // Now cut the cylinder. Keep in mind that the data provided always uses a cylinder with xy planes and z as
        // height.
        auto a = bb.first;
        auto b = bb.second;
        // How much?
        auto frac = 1.0 - (cutPerc / 100.0);
        auto fracZ = 1.0 - (cutPercZ / 100.0);

        // Distance along X
        auto dX = b[0] - a[0];
        // Distance along y
        auto dY = b[1] - a[1];
        // Distance along z
        auto dZ = b[2] - a[2];

        // Define the cutBB for the network data
        return std::make_pair(Vec3({{a[0] + (0.5 * frac * dX), a[1] + (0.5 * frac * dY), a[2] + (0.5 * fracZ * dZ)}}),
                              Vec3({{b[0] - (0.5 * frac * dX), b[1] - (0.5 * frac * dY), b[2] - (0.5 * fracZ * dZ)}}));
    }

    /**
     * Calculate the elliptic cylinder volume given by a BB
     *
     * \param bb the bb
     *
     * \return the volume
     */
    auto cylinderVolume(const BoundingBox& bb)
    {
        // Ellipse Cylinder Volume: pi * a * b * h

        // The volume of the whole cylinder (uncut)
        auto a = 0.5 * (bb.second[0] - bb.first[0]);
        auto b = 0.5 * (bb.second[1] - bb.first[1]);
        auto h = bb.second[2] - bb.first[2];
        return pi * a * b * h;
    }

    /**
     * Calculate the volume given by a BB
     *
     * \param bb the bb
     * \param cylindric true if the bb represents a cylinder
     *
     * \return the volume
     */
    auto volume(const BoundingBox& bb, bool cylindric)
    {
        if (cylindric)
        {
            return cylinderVolume(bb);
        }

        return (bb.second[0] - bb.first[0]) * (bb.second[1] - bb.first[1]) * (bb.second[2] - bb.first[2]);
    }

    /**
     * Check the normalized distance to a cylinder defined by a bounding box.
     *
     * \param p the point to check
     * \param bb the BB representing the bounding box with h in z-direction
     *
     * \return normalized distance to center as pair. First is radial distance, second is z-distance to center.
     */
    auto pointCylinderDistanceNorm(const Vec3& p, const BoundingBox& bb)
    {
        // is inside the cylinder of cutBB?
        auto a = bb.first;
        auto b = bb.second;

        // Get the centerline of the cylinder
        auto rX = 0.5 * (b[0] - a[0]);
        auto rY = 0.5 * (b[1] - a[1]);
        auto rZ = 0.5 * (b[2] - a[2]);
        auto mX = a[0] + rX;
        auto mY = a[1] + rY;
        auto mZ = a[2] + rZ;

        // Normalize the point to a unit cylinder
        auto npX = (p[0] - mX) / rX;
        auto npY = (p[1] - mY) / rY;
        auto npZ = (p[2] - mZ) / rZ;

        // Get the radius that would present a circle in m through npX,npY
        auto r = std::sqrt(std::pow(npX, 2.0) + std::pow(npY, 2.0));

        return std::make_pair(r, std::abs(npZ));
    }

    /**
     * Calc an intersection with the cylinder between p1 and p2. P2 needs to be outside. The intersection is done via
     * approximation, not exact.
     *
     * \param p1 inside point
     * \param p2 outside point
     * \param count the recursion depth.
     * \param bb the bounding box representing a cylinder
     *
     * \return the point on the cylinder surface
     */
    Vec3 cylinderIntersection(const Vec3& p1, const Vec3& p2, const BoundingBox& bb, size_t count = 0)
    {
        // Get center
        auto center = add(p2, scale(difference(p1, p2), 0.5));

        // Since the lines where beautified using splines, they are  sampled densely. The center point is probably
        // almost near the border.
        auto d = pointCylinderDistanceNorm(center, bb);
        auto off = 1.0 - d.first;
        if ((off < 0.001) || count >= 3)
        {
            return center;
        }

        // Probably find a better approx:
        if (d.first < 1.0)
        {
            return cylinderIntersection(center, p2, bb, count + 1);
        }
        else if (d.first == 1.0)
        {
            return center;
        }
        else if (d.first > 1.0)
        {
            return cylinderIntersection(p1, center, bb, count + 1);
        }

        return center;
    }

    /**
     * Derive a network from plain, indexed line data.
     *
     * \param points the points to use
     * \param lines the lines to use
     * \param radii of each line
     * \param splitAtDiffRadii if true, different radii cause a segment to be split.
     * \param mergeCloseBranchPoints if true, the algorithm merges branchpoints not more far away than maxMergeDistance
     * \param maxMergeDistance branchpoints closer than this distance will be merged if mergeCloseBranchPoints is true
     *
     * \return a vessel segments structure.
     */
    VesselSegments deriveNetwork(const decltype(VesselNetwork::points)& points,
                                 const decltype(VesselNetwork::lines)& lines,
                                 const decltype(VesselNetwork::radii)& radii, bool splitAtDiffRadii,
                                 bool mergeCloseBranchPoints, Real maxMergeDistance)
    {
        VesselSegments vessels;
        vessels.points = points;

        // We need some more information first:

        // Do as long as we have non-assigned lines. We mark each assigned line as "assigned" and skip them.
        std::vector< bool > lineAssigned(lines.size(), false);

        LogD << "Calculating degree of each point." << LogEnd;

        // Gather information about each point's degree - how many lines are adjacent to this point? Index matches point
        // index.
        vessels.degrees.resize(vessels.points.size(), 0);

        // Btw: this is our chance to build an inverse points-to-lines map (adjacency list):
        std::vector< std::vector< VesselNetwork::IndexType > > linesByPoint(vessels.points.size());

        std::vector< std::set< Real > > pointRadii;
        pointRadii.resize(vessels.points.size());

        vessels.volumeCapillaries = 0.0;
        vessels.volumeNonCapillaries = 0.0;

        // Iterate each line:
        for (auto lineIdx : range(lines.size()))
        {
            auto p1Idx = lines[lineIdx].first;
            auto p1 = points[p1Idx];
            auto p2Idx = lines[lineIdx].second;
            auto p2 = points[p2Idx];

            // NOTE: keep in mind that this will be equal to linesByPoint[ p1( or 2)Idx ].size(). But as we might need
            // this value quite often, we store it separately to avoid multiple size() calls.
            vessels.degrees[p1Idx]++;
            vessels.degrees[p2Idx]++;

            // Build adjacency list
            linesByPoint[p1Idx].push_back(lineIdx);
            linesByPoint[p2Idx].push_back(lineIdx);

            // Build adjacent radii list
            pointRadii[p1Idx].insert(radii[lineIdx]);
            pointRadii[p2Idx].insert(radii[lineIdx]);

            // Also sum up volumes of the lines
            auto dist = distance(p1, p2);
            auto vol = pi * std::pow(radii[lineIdx], 2.0) * dist;
            vessels.volumeCapillaries += isCapillary(radii[lineIdx]) ? vol : 0.0;
            vessels.volumeNonCapillaries += !isCapillary(radii[lineIdx]) ? vol : 0.0;
        }

        LogD << "Calculating the branch points." << LogEnd;

        // A criterion is required to check whether a point is assumed to be an end-point of a segment.
        auto isEndPoint = [&](const VesselNetwork::IndexType& pointIdx) {
            // A point is NOT an endpoint if it is adjacent to 2 lines and both lines have the same radius. ->
            // accordingly, it is an endpoint if this is not true.
            auto degreeCriterion = (vessels.degrees[pointIdx] != 2);
            auto radiusCriterion = (pointRadii[pointIdx].size() != 1);
            return degreeCriterion || (splitAtDiffRadii && radiusCriterion);
        };

        // Lets use an in-line function for this. If we would write our own function, it would need a whole bunch of
        // additional parameters for the lines vector, points vector linesByPoint .... and so on.
        //
        // This method also associates the originally given radius to each segment. It compares the radius to the
        // expected radius, as the data SHOULD provide a single radius for all lines of a single segment. Fails if this
        // is not the case!
        //
        // The function returns true if there has been an issue with multiple radii at a single segment.
        std::function< Real(bool, VesselSegments::Segment&, const VesselNetwork::IndexType&,
                            const VesselNetwork::IndexType&, VesselSegments::SegmentRadii&) >
            traceSegment =
                // NOTE: avoid auto declarations here as they remove the reference
            [&](bool front, VesselSegments::Segment& segment, const VesselNetwork::IndexType& previousPointIdx,
                const VesselNetwork::IndexType& lineIdx, VesselSegments::SegmentRadii& segRadii) {
                // Mark as done
                lineAssigned[lineIdx] = true;

                // Get points
                auto p1Idx = lines[lineIdx].first;
                auto p2Idx = lines[lineIdx].second;
                auto dist = distance(vessels.points[p1Idx], vessels.points[p2Idx]);

                // Which point is the next to add?
                auto nextPointIdx = (p1Idx == previousPointIdx) ? p2Idx : p1Idx;

                // Some datasets contain cycles. To prevent endless looping, check if the next point was visited
                // already.
                bool loopFound = (std::find(segment.begin(), segment.end(), nextPointIdx) != segment.end());

                // There are three interesting cases now. The next point is
                // 1) a branch point ( degree > 2 )
                // 2) an end point ( degree == 1 )
                // 3) an intermediate point between the current line and the next line ( degree == 2 )

                // All cases cause the current point to be added:
                segRadii[radii[lineIdx]]++;
                if (front)
                {
                    segment.push_front(nextPointIdx);
                }
                else
                {
                    segment.push_back(nextPointIdx);
                }

                // We check AFTER adding the point. This way, the first and last point of the segment are equal (closed
                // loop). This is more accurate for calculating the segment length.
                if (loopFound)
                {
                    LogW << "Loop found. Gracefully aborting trace." << LogEnd;
                    return 0.0;
                }

                // In cases 1 and 2, recursion stops. So handle case 3, stop in all other cases.
                if (!isEndPoint(nextPointIdx)) // intermediate point?
                {
                    // Get next line index:
                    auto nextLineIdx = (linesByPoint[nextPointIdx].front() == lineIdx)
                                           ? linesByPoint[nextPointIdx].back()
                                           : linesByPoint[nextPointIdx].front();

                    // Recursive call ...
                    return dist + traceSegment(front, segment, nextPointIdx, nextLineIdx, segRadii);
                }

                // Done. We reached an endpoint.
                return dist;
            };

        // Keep in mind that this routine might change the order of points in a line. So a line defined by (a,b), might
        // be represented here as (b,a) as we cannot make any assumptions on the ordering of points in the source data.

        // Iterate. Trace each unassigned line until all lines have been traces/assigned.
        for (auto lineIdx : range(lines.size()))
        {
            // If the line was assigned somewhere, skip tracing it.
            if (lineAssigned[lineIdx])
            {
                continue;
            }

            // Create a new segment
            VesselSegments::Segment segment;

            // Trace in both directions:
            VesselSegments::SegmentRadii segRadii;
            auto l = traceSegment(true, segment, lines[lineIdx].second, lineIdx, segRadii);
            l += traceSegment(false, segment, lines[lineIdx].first, lineIdx, segRadii);

            // Add to list
            vessels.segments.push_back(segment);
            vessels.removed.push_back(false);

            // Store radii too
            vessels.radii.push_back(segRadii);

            // and length
            vessels.length.push_back(l);
        }

        // Clean up some special cases of parallel segments with only two points. There are a lot of cases, where
        // several lines connect the same points. Remove them.
        size_t remCount = 0;
        for (auto segmentIdx : range(vessels.segments.size()))
        {
            if (vessels.removed[segmentIdx])
            {
                continue;
            }
            if (vessels.segments[segmentIdx].size() == 2)
            {
                for (auto segmentIdx2 : range(vessels.segments.size()))
                {
                    // skip self
                    if (segmentIdx2 == segmentIdx)
                    {
                        continue;
                    }
                    if (vessels.removed[segmentIdx2])
                    {
                        continue;
                    }
                    if (vessels.segments[segmentIdx2].size() == 2)
                    {
                        auto c1 = (vessels.segments[segmentIdx].front() == vessels.segments[segmentIdx2].back());
                        auto c2 = (vessels.segments[segmentIdx].back() == vessels.segments[segmentIdx2].front());
                        auto same = c1 && c2;

                        auto c3 = (vessels.segments[segmentIdx].front() == vessels.segments[segmentIdx2].front());
                        auto c4 = (vessels.segments[segmentIdx].back() == vessels.segments[segmentIdx2].back());
                        same = same || (c3 && c4);

                        vessels.removed[segmentIdx2] = same;
                        remCount += same ? 1 : 0;
                    }
                }
            }
        }

        LogI << "Found " << vessels.segments.size() << " vessel segments. Removed " << remCount
             << " as they where equal." << LogEnd;

        if (mergeCloseBranchPoints)
        {
            LogI << "Merging short segments. Max distance = " << maxMergeDistance << LogEnd;

            // Check each segment:
            size_t mergeCount = 0;
            size_t numMergeIterations = 0;
            size_t numRemoved = 0;
            bool wasMerging = false;
            do
            {
                wasMerging = false;
                numMergeIterations++;
                for (auto segmentIdx : range(vessels.segments.size()))
                {
                    // Ignore segments that have been removed already
                    if (vessels.removed[segmentIdx])
                    {
                        continue;
                    }

                    auto l = vessels.length[segmentIdx];
                    auto frontIdx = vessels.segments[segmentIdx].front();
                    auto backIdx = vessels.segments[segmentIdx].back();

                    // length too short and between two branch points?
                    if ((l < maxMergeDistance) && (vessels.degrees[frontIdx] > 2) && (vessels.degrees[backIdx] > 2))
                    {
                        // Merge.
                        mergeCount++;
                        wasMerging = true;
                        vessels.removed[segmentIdx] = true;

                        // 1) get new point coordinate:
                        auto newPoint = midpoint(vessels.points[frontIdx], vessels.points[backIdx]);

                        // 2) add new point to points list
                        auto newIdx = vessels.points.size();
                        vessels.points.push_back(newPoint);
                        vessels.degrees.push_back(vessels.degrees[frontIdx] + vessels.degrees[backIdx] - 2);
                        vessels.degrees[frontIdx] = 0;
                        vessels.degrees[backIdx] = 0;

                        // 3) tell all involved segments that there is a new point
                        for (auto segmentIdxInvolved : range(vessels.segments.size()))
                        {
                            // skip self
                            if (segmentIdxInvolved == segmentIdx)
                            {
                                continue;
                            }

                            // skip removed segments
                            if (vessels.removed[segmentIdxInvolved])
                            {
                                continue;
                            }

                            auto frontIdxInvolved = vessels.segments[segmentIdxInvolved].front();
                            auto backIdxInvolved = vessels.segments[segmentIdxInvolved].back();

                            // Change all the involved segments and update their lengths
                            Real oldD = 0.0;
                            Real newD = 0.0;
                            size_t numEndsChanged = 0;
                            bool needLengthUpdate = false;
                            if ((frontIdxInvolved == frontIdx) || (frontIdxInvolved == backIdx))
                            {
                                // The new length between the old front and its neighbour element
                                auto p1Idx = frontIdxInvolved;
                                auto p2Idx = *(std::next(vessels.segments[segmentIdxInvolved].begin()));
                                oldD += distance(vessels.points[p1Idx], vessels.points[p2Idx]);
                                newD += distance(vessels.points[newIdx], vessels.points[p2Idx]);

                                vessels.segments[segmentIdxInvolved].front() = newIdx;
                                needLengthUpdate = true;
                                numEndsChanged++;
                            }
                            if ((backIdxInvolved == frontIdx) || (backIdxInvolved == backIdx))
                            {
                                // The new length between the old front and its neighbour element
                                auto p1Idx = backIdxInvolved;
                                auto p2Idx = *(std::prev(vessels.segments[segmentIdxInvolved].end(), 2));
                                oldD += distance(vessels.points[p1Idx], vessels.points[p2Idx]);
                                newD += distance(vessels.points[newIdx], vessels.points[p2Idx]);

                                vessels.segments[segmentIdxInvolved].back() = newIdx;
                                needLengthUpdate = true;
                                numEndsChanged++;
                            }

                            // There is a special case to handle: Some segments might now get obsolete, as they now use
                            // the newIdx point as start and end. Remove them.
                            if ((vessels.segments[segmentIdxInvolved].size() == 2) && (numEndsChanged == 2))
                            {
                                vessels.removed[segmentIdxInvolved] = true;
                                numRemoved++;
                            }
                            else if (needLengthUpdate)
                            {
                                auto currentLength = vessels.length[segmentIdxInvolved];
                                vessels.length[segmentIdxInvolved] = currentLength - oldD + newD;
                                if (vessels.length[segmentIdxInvolved] < 0.0)
                                {
                                    LogE << vessels.segments[segmentIdxInvolved].size() << " - " << frontIdxInvolved
                                         << " --- " << backIdxInvolved
                                         << " has negative length! This should not happen. Current = " << currentLength
                                         << " minus old " << oldD << " plus new " << newD << LogEnd;
                                }
                            }
                        }
                    }
                }
                // Thats it ...
            } while (wasMerging);

            LogD << "Merged " << mergeCount << " segments in " << numMergeIterations - 1 << " iterations. Removed "
                 << numRemoved << " due to equality to merged vessels." << LogEnd;
        }

        LogI << "Building segment adjacency list." << LogEnd;
        for (auto segmentIdx : range(vessels.segments.size()))
        {
            if (vessels.removed[segmentIdx])
            {
                continue;
            }
            for (auto segmentIdx2 : range(vessels.segments.size()))
            {
                // skip self
                if (segmentIdx2 == segmentIdx)
                {
                    continue;
                }
                if (vessels.removed[segmentIdx2])
                {
                    continue;
                }

                // They are adjacent of start/end indices match somehow
                auto c1 = (vessels.segments[segmentIdx].front() == vessels.segments[segmentIdx2].back());
                auto c3 = (vessels.segments[segmentIdx].front() == vessels.segments[segmentIdx2].front());
                if (c1 || c3)
                {
                    vessels.adjacency[vessels.segments[segmentIdx].front()].insert(segmentIdx2);
                }

                auto c2 = (vessels.segments[segmentIdx].back() == vessels.segments[segmentIdx2].front());
                auto c4 = (vessels.segments[segmentIdx].back() == vessels.segments[segmentIdx2].back());
                if (c2 || c4)
                {
                    vessels.adjacency[vessels.segments[segmentIdx].back()].insert(segmentIdx2);
                }
            }
        }

        // Calculate new degrees since multiple merges might confuse the upper degree calculation
        for (auto item : vessels.adjacency)
        {
            auto pIdx = item.first;
            auto&& adj = item.second;

            auto diff = static_cast< int64_t >(adj.size()) - static_cast< int64_t >(vessels.degrees[pIdx]);
            if (diff != 0)
            {
                // LogW << pIdx << " has different degrees in adjacency - calculated: " << adj.size() << " - " <<
                // vessels.degrees[ pIdx ] << " = "
                //      << diff << LogEnd;
                //  vessels.degrees[ pIdx ] = adj.size();
            }
        }

        return vessels;
    }

    /**
     * Check if a point is inside a given volume
     *
     * \param p the point to check
     * \param cylindric if true, assume a cylinder instead of a cube in bb.
     *
     * \return true if inside.
     */
    bool pointInVolume(const Vec3& p, const BoundingBox& bb, bool cylindric)
    {
        if (cylindric)
        {
            // Check if this new cricle is at most as large as the desired one
            auto d = pointCylinderDistanceNorm(p, bb);
            return ((d.first <= 1.0) && (d.second <= 1.0));
        }

        // Upper inclusive
        return (bb.first <= p) && (p <= bb.second);
    }

    /**
     * Cut the given volume and mask everything outside the given bounding volume.
     *
     * \param srcVolume volume to cut
     * \param bb target BB representing a volume with height = z
     * \param cylindric if true, assume a cylinder instead of a cube in bb.
     * \param outsideValue the value to write to voxels outside the volume. Default is ValueType()
     *
     * \return volume and the count of voxels inside the volume.
     */
    template < typename ValueType >
    std::pair< Volume< ValueType >, size_t > cutVolume(const Volume< ValueType >& srcVolume, const BoundingBox& bb,
                                                       bool cylindric, ValueType outsideValue = ValueType())
    {
        LogI << "Cutting volume." << LogEnd;
        Volume< ValueType > volume(bb.second[0] - bb.first[0], bb.second[1] - bb.first[1], bb.second[2] - bb.first[2]);

        // Just create a plain copy. Be aware, that this does not include the cylinder test. This needs to be done
        // during calculations.
        size_t inVoxelCount = 0;
        for (auto z : range(bb.first[2], bb.second[2]))
        {
            for (auto y : range(bb.first[1], bb.second[1]))
            {
                for (auto x : range(bb.first[0], bb.second[0]))
                {
                    auto srcIdx = srcVolume.index(x, y, z);
                    auto targetIdx = volume.index(x - bb.first[0], y - bb.first[1], z - bb.first[2]);

                    // Valid in cylinder?
                    auto p = Vec3({{static_cast< Real >(x), static_cast< Real >(y), static_cast< Real >(z)}});
                    int isIn = pointInVolume(p, bb, cylindric) ? 1 : 0;
                    inVoxelCount += isIn;

                    volume.data[targetIdx] = isIn ? srcVolume.data[srcIdx] : outsideValue;
                }
            }
        }

        return std::make_pair(volume, inVoxelCount);
    }

    /**
     * Process each voxel.
     *
     * \tparam FunctorType the functor to use per voxel
     * \param bb the bounding box describing the cylinder or box
     * \param cylindric if true, assume a cylinder instead of a cube in bb.
     * \param functor the functor instance to call
     *
     * \return the number of processed voxels
     */
    template < typename FunctorType >
    size_t processVoxels(const BoundingBox& bb, bool cylindric, FunctorType functor)
    {
        size_t usedVoxels = 0;
        for (auto z : range(bb.first[2], bb.second[2]))
        {
            for (auto y : range(bb.first[1], bb.second[1]))
            {
                for (auto x : range(bb.first[0], bb.second[0]))
                {
                    auto p = Vec3({{static_cast< Real >(x), static_cast< Real >(y), static_cast< Real >(z)}});
                    if (pointInVolume(p, bb, cylindric))
                    {
                        functor(x, y, z);
                        usedVoxels++;
                    }
                }
            }
        }
        return usedVoxels;
    }

    /**
     * Calculate a discrete bin number of a given value.
     *
     * @tparam ValueType the type used for values
     * @tparam SizeType the size type used to define bin sizes
     * \param value the value to get the bin number of
     * \param numBins the number of bins
     * \param minVal the smallest value allowed
     * \param maxVal the largest value allowed
     *
     * \return the bin number
     */
    template < typename ValueType, typename SizeType >
    int calcBinNum(const ValueType& value, const SizeType& numBins, const ValueType& minVal, const ValueType& maxVal)
    {
        return static_cast< int >(std::floor(static_cast< Real >(numBins) * (value - minVal) / (maxVal - minVal)));
    }

    /**
     * Do actual analysis of the data. If something fails, an exception will be thrown.
     *
     * \param unorderedVesselData the vessel network data. It is assumed that this data is unordered.
     * \param maskVolume the vessel volume data.
     * \param dtVolume the distance transformed volume. Use nullptr to skip EVD calculations
     * \param outputDir the directory where to write everything. Ensure a trailing / or \!
     * \param settings the settings to use
     */
    void analyze(std::string outputDir, const VesselNetwork& unorderedVesselData, ConstSPtr< MaskVolume > maskVolume,
                 ConstSPtr< DTVolume > dtVolume, const Settings& settings)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 1a - Given
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // This will keep the data we calculate. Will be written to file eventually.
        SegmentData localSegmentData;
        BranchPointData localBranchPointData;
        VolumeData globalVolumeData;
        EVDData globalEVDData;

        // As the typical programmer is lazy as sh**, we provide some easy access to the members of the vessel network
        const decltype(unorderedVesselData.points)& origPoints = unorderedVesselData.points;
        const decltype(unorderedVesselData.lines)& origLines = unorderedVesselData.lines;
        const decltype(unorderedVesselData.radii)& origRadii = unorderedVesselData.radii;
        const BB& bb = unorderedVesselData.boundingBox;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 1b - Assumptions
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        bool isCylindric = settings.cylindric;

        // According to Thomas, the data needs to be cut by 50% along the cylinder radius
        Real networkCutXY = settings.networkCutXY;
        Real networkCutZ = settings.networkCutZ;

        // Similar for volume data:
        Real volumeCutXYVVF = settings.volumeCutXYVVF;
        Real volumeCutZVVF = settings.volumeCutZVVF;

        // Similar for "locally influenced" volume data:
        // IF you change this to be something different from the network cut, also adapt the below calculations that use
        // the network volume
        Real volumeCutXYVVFLocalized = networkCutXY;
        Real volumeCutZVVFLocalized = networkCutZ;

        // Similar for evd volume data:
        Real volumeCutXYEVD = settings.volumeCutXYEVD;
        Real volumeCutZEVD = settings.volumeCutZEVD;

        // The original publication ignored branchpoints in the statistics with a certain degree. Turn this on/off here.
        bool restrictDegree = settings.restrictDegree;
        bool mapLargerToMax = settings.mapLargerToMax;
        size_t maxBranchpointDegree = settings.maxBranchpointDegree;

        // Do a merging step to merge branch points that are very close?
        bool mergeCloseBranchPoints = settings.mergeCloseBranchPoints;
        Real maxMergeDistance = settings.maxMergeDistance;

        // Do we split segments in case of multiple radii per segement?
        bool splitSegmentsAtDiffRadii = false;

        // We store the diameter per segment length in a discrete length-span: Define it here:
        Real diameterPerLengthMin = settings.diameterPerLengthMin;
        Real diameterPerLengthMax = settings.diameterPerLengthMax;
        auto diameterPerLengthNumBins = settings.diameterPerLengthNumBins;

        // Do we define a minimal EVD value?
        Real minEVD = settings.minEVD;
        bool useMinEVD = settings.useMinEVD;
        Real maxEVD = settings.maxEVD;
        bool useMaxEVD = settings.useMaxEVD;

        // Define a discretization scheme for EVD. This defines the histogram used to categorize EVD values.
        Real evdHistMin = settings.evdHistMin;
        Real evdHistMax = settings.evdHistMax;
        auto evdNumBins = settings.evdNumBins;

        // Scale the EVD. The Evd was calculated in voxel space. Use the voxel size in the volume as scaler here:
        auto evdScale = settings.pixelSize;
        // Real evdScale = 0.725586;
        // Real evdScale = 0.729492;

        bool useVVF = (maskVolume != nullptr);
        bool useEVD = (dtVolume != nullptr);

        // We need to apply the cut first. On the cut data, we calculate the statistics
        LogD << "Bounding Box: [ ( " << bb.first[0] << ", " << bb.first[1] << ", " << bb.first[2] << " ), ( "
             << bb.second[0] << ", " << bb.second[1] << ", " << bb.second[2] << " ) ]" << LogEnd;

        // We need to apply the cut first. On the cut data, we calculate the statistics
        auto cutNetworkBB = cutBoundingBox(bb, networkCutXY, networkCutZ);
        LogD << "Cut Network Bounding Box: [ ( " << cutNetworkBB.first[0] << ", " << cutNetworkBB.first[1] << ", "
             << cutNetworkBB.first[2] << " ), ( " << cutNetworkBB.second[0] << ", " << cutNetworkBB.second[1] << ", "
             << cutNetworkBB.second[2] << " ) ]" << LogEnd;

        // A criterium defining whether a point is inside our target domain
        auto pointOK = [&](auto&& p) { return pointInVolume(p, cutNetworkBB, isCylindric); };

        // Important: this code avoids actually preprocessing the line data. Instead, the below code uses the
        // ignore-predicate to check at each case whether to include the point into the statistics or not.

        // The whole volume. The volume is in µm^3. We require mm^3 though.
        Real wholeDomainVolume = volume(bb, isCylindric) / std::pow(1000.0, 3);

        // The cut volume. The volume is in µm^3. We require mm^3 though.
        Real networkCutDomainVolumeMicroMeterCube = volume(cutNetworkBB, isCylindric);
        Real networkCutDomainVolume = networkCutDomainVolumeMicroMeterCube / std::pow(1000.0, 3);

        LogD << "Original volume of " << wholeDomainVolume << " mm^3" << LogEnd;
        LogD << "Cut volume of " << networkCutDomainVolume << " mm^3" << LogEnd;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 1c - cutting the data
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        LogI << "Cutting network." << LogEnd;

        decltype(unorderedVesselData.points) points;
        decltype(unorderedVesselData.lines) lines;
        decltype(unorderedVesselData.radii) radii;

        // keep track of the new index of a point.
        std::vector< int64_t > newPointIndex(origPoints.size(), -1);

        // For each line
        for (auto lineIdx : range(origLines.size()))
        {
            // Copy the points if not done already:
            auto p1Idx = origLines[lineIdx].first;
            auto p2Idx = origLines[lineIdx].second;

            auto p1 = origPoints[p1Idx];
            auto p2 = origPoints[p2Idx];

            auto p1NewIdx = newPointIndex[p1Idx];
            auto p2NewIdx = newPointIndex[p2Idx];

            // Handle three cases:
            if (!settings.cylindric || settings.preCroppedNetwork || (pointOK(p1) && pointOK(p2)))
            {
                // Case 1: all points inside the volume
                // handle?
            }
            else if (pointOK(p1) && !pointOK(p2))
            {
                // Case 2: one point outside

                // Interpolate the other one
                p2 = cylinderIntersection(p1, p2, cutNetworkBB);
                // This is definitely a new point.
                p2NewIdx = -1;
            }
            else if (pointOK(p2) && !pointOK(p1))
            {
                // Case 2: one point outside, the other way around

                // Interpolate the other one
                p1 = cylinderIntersection(p2, p1, cutNetworkBB);
                // This is definitely a new point.
                p1NewIdx = -1;
            }
            else
            {
                // Case 3: both points are outside the new volume. Perfect. Nothing to do .. Ignore the whole line and
                // its points. Ignore the fact that a long line could cross the volume. The volume is much larger than
                // the longest segments.
                continue;
            }

            // Add the points of not yet existing
            if (p1NewIdx == -1)
            {
                p1NewIdx = points.size();
                points.push_back(p1);
                newPointIndex[p1Idx] = p1NewIdx;
            }
            if (p2NewIdx == -1)
            {
                p2NewIdx = points.size();
                points.push_back(p2);
                newPointIndex[p2Idx] = p2NewIdx;
            }

            // Add the line
            lines.push_back(std::make_pair(p1NewIdx, p2NewIdx));
            radii.push_back(origRadii[lineIdx]);
        }

        LogD << "Original Network Info: " << origPoints.size() << " points in " << origLines.size() << " lines."
             << LogEnd;
        LogD << "Cut Network Info: " << points.size() << " points in " << lines.size() << " lines." << LogEnd;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 2 - Derive the network.
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        LogI << "Reconstructing vessel segments." << LogEnd;
        // And this one contains the derived real network.
        auto vessels =
            deriveNetwork(points, lines, radii, splitSegmentsAtDiffRadii, mergeCloseBranchPoints, maxMergeDistance);

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 3 - Analyze network and derive segment data
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Count.
        size_t numSegments = 0;
        size_t numSegmentsNotBetweenBranchpoints = 0;
        size_t numCapillarySegments = 0;
        size_t numNonCapillarySegments = 0;

        // Collect the volumes
        Real volumeCap = 0.0;
        Real volumeNonCap = 0.0;

        // Binned segment lenghts
        Real diameterPerLengthBinWidth =
            (diameterPerLengthMax - diameterPerLengthMin) / static_cast< Real >(diameterPerLengthNumBins);
        std::vector< std::vector< Real > > diameterPerLenght(diameterPerLengthNumBins);

        localSegmentData.anglesAdjacency.reserve(vessels.segments.size() * 2);

        struct AdjacentBifVector
        {
            Vec3 gradient;
            Real radius;
        };
        std::map< size_t, std::vector< AdjacentBifVector > > adjacentBifVectors;

        // auto vecToMajorDirection = []( const auto& a )
        // {
        //     int dir = 0;
        //     if( ( a[1] > a[0] ) && ( a[1] > a[2] ) )
        //     {
        //         dir = 1;
        //     }
        //     else if( ( a[2] > a[0] ) && ( a[2] > a[1] ) )
        //     {
        //         dir = 2;
        //     }
        //
        //     return dir;
        // };

        // We have a network -> copy its data and calc if needed
        for (auto segmentIdx : range(vessels.segments.size()))
        {
            if (vessels.removed[segmentIdx])
            {
                continue;
            }

            auto&& segment = vessels.segments[segmentIdx];

            if ((vessels.degrees[segment.front()] < 3) || (vessels.degrees[segment.back()] < 3))
            {
                numSegmentsNotBetweenBranchpoints++;
                // continue;
            }

            // All containers have the same size. As we push back, but skip some segments, segmentIDX and the real index
            // in these containers may be different -> store in a map
            auto realIdx = localSegmentData.directLengths.size();
            localSegmentData.idxMap[segmentIdx] = realIdx;

            // Length
            auto realLength = vessels.length[segmentIdx];
            localSegmentData.lengths.push_back(realLength);

            // Direct distance
            auto directDistance = distance(vessels.points[segment.front()], vessels.points[segment.back()]);
            localSegmentData.directLengths.push_back(directDistance);

            // Tortuosity and direction
            if (directDistance > 0.0)
            {
                auto tortuosity = realLength / directDistance;
                if (segment.front() == segment.back())
                {
                    LogI << "Loop found. Tortuosity is not written for this segment." << LogEnd;
                }
                else
                {
                    localSegmentData.tortuosity.push_back(tortuosity);
                }
            }
            else
            {
                LogI << "Loop found. Tortuosity is not written for this segment." << LogEnd;
            }

            // Radius and Diameter ... well this can be a problem. How to handle multiple radii for a single segment?
            size_t sumRadii = 0;
            auto radius = std::accumulate(vessels.radii[segmentIdx].begin(), vessels.radii[segmentIdx].end(), 0.0,
                                          [&](auto&& a, auto&& b) {
                                              sumRadii += b.second;
                                              return a + b.first * static_cast< Real >(b.second);
                                          }) /
                          static_cast< Real >(sumRadii);

            auto isACapillary = isCapillary(radius);

            // v = pi * r^2 * h
            Real vol = pi * std::pow(radius, 2.0) * realLength;

            localSegmentData.radiiSingle[segmentIdx] = radius;
            localSegmentData.radius.push_back(radius);
            localSegmentData.diameter.push_back(radius * 2.0);
            localSegmentData.isCapillary.push_back(isACapillary);
            localSegmentData.volume.push_back(vol);

            if ((directDistance > 0.0) && (segment.front() != segment.back()))
            {
                // Orientation vector for segments is the vector between end and start
                auto segGradient = difference(vessels.points[segment.back()], vessels.points[segment.front()]);
                auto segLen = length(segGradient);
                auto segNormalizedGradient = scale(segGradient, 1.0 / Real(segLen));

                if (segLen > 0.0)
                {
                    // segGradient is from front (p1) to back (p2)
                    if (vessels.degrees[segment.front()] > 2)
                    {
                        adjacentBifVectors[segment.front()].emplace_back(
                            AdjacentBifVector{segNormalizedGradient, radius});
                    }
                    if (vessels.degrees[segment.back()] > 2)
                    {
                        adjacentBifVectors[segment.back()].emplace_back(
                            AdjacentBifVector{scale(segNormalizedGradient, -1.0), radius});
                    }
                }
            }

            // Some stats
            numCapillarySegments += isACapillary ? 1 : 0;
            numNonCapillarySegments += !isACapillary ? 1 : 0;
            numSegments++;

            // Sum up the volumes
            volumeCap += isACapillary ? vol : 0.0;
            volumeNonCap += !isACapillary ? vol : 0.0;

            // It might be of interest to get the mean diameter for a segment length. So, collect the diameters
            auto binNum = calcBinNum(realLength, diameterPerLengthNumBins, diameterPerLengthMin, diameterPerLengthMax);
            if ((binNum < 0) || (binNum >= diameterPerLengthNumBins))
            {
                continue;
            }
            diameterPerLenght[binNum].push_back(radius * 2.0);
        }

        localSegmentData.segmentDensity = static_cast< Real >(numSegments) / networkCutDomainVolume;
        localSegmentData.segmentDensityCapillary = static_cast< Real >(numCapillarySegments) / networkCutDomainVolume;
        localSegmentData.segmentDensityNonCapillary =
            static_cast< Real >(numNonCapillarySegments) / networkCutDomainVolume;

        localSegmentData.numSegments = numSegments;
        localSegmentData.numCapillarySegments = numCapillarySegments;
        localSegmentData.numNonCapillarySegments = numNonCapillarySegments;

        // Accumulate bins and calc mean per bin
        localSegmentData.meanDiameterPerLength.resize(diameterPerLengthNumBins);
        {
            size_t binNum = 0;
            std::generate(
                localSegmentData.meanDiameterPerLength.begin(), localSegmentData.meanDiameterPerLength.end(), [&]() {
                    // Calc mean:
                    auto mean =
                        std::accumulate(diameterPerLenght[binNum].begin(), diameterPerLenght[binNum].end(), Real());
                    if (diameterPerLenght[binNum].size() == 0)
                    {
                        return 0.0;
                    }
                    mean /= static_cast< Real >(diameterPerLenght[binNum].size());
                    binNum++;

                    return mean;
                });
        }

        localSegmentData.meanDiameterPerLengthBinCenters.resize(diameterPerLengthNumBins);
        {
            size_t binNum = 0;
            std::generate(localSegmentData.meanDiameterPerLengthBinCenters.begin(),
                          localSegmentData.meanDiameterPerLengthBinCenters.end(), [&]() {
                              Real center = diameterPerLengthMin + (diameterPerLengthBinWidth / 2.0) +
                                            (static_cast< Real >(binNum) * diameterPerLengthBinWidth);
                              binNum++;
                              return center;
                          });
        }

        LogD << numSegmentsNotBetweenBranchpoints << " of " << numSegments << " segments are not between branchpoints."
             << LogEnd;

        // Calculate segment angles at branch points
        for (auto item : adjacentBifVectors)
        {
            // const auto& pIdx = item.first;
            auto& adjBifVectors = item.second;

            std::sort(adjBifVectors.begin(), adjBifVectors.end(),
                      [](const auto& a, const auto& b) { return a.radius > b.radius; });

            for (auto adjVecIdx : range(adjBifVectors.size()))
            {
                const auto& refVec = adjBifVectors[adjVecIdx];
                for (auto relVecIdx : range(adjVecIdx + 1, adjBifVectors.size()))
                {
                    const auto& relVec = adjBifVectors[relVecIdx];
                    auto angle = degree(refVec.gradient, relVec.gradient);

                    localSegmentData.anglesAdjacency.addByDiameter(relVec.radius * 2.0, angle,
                                                                   isCapillary(relVec.radius));
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 4 - Analyze network and derive branchpoint data
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        auto highestDegree = *std::max_element(vessels.degrees.begin(), vessels.degrees.end());

        // The highest degree to include into the statistics
        auto maxDegree = restrictDegree ? std::min(highestDegree, maxBranchpointDegree) : highestDegree;

        // Reserve some space -> reduce allocations during loop
        localBranchPointData.degrees.reserve(vessels.points.size());

        // This collects some info used for user debug
        std::map< size_t, size_t > branchPointCountByDegree;

        // As points are shared, we need to keep track which point was handled already.
        std::map< size_t, bool > handledPoint;

        // We have a network -> each segments start/end is either a branch point or a "finger" - endpoint
        for (auto segmentIdx : range(vessels.segments.size()))
        {
            if (vessels.removed[segmentIdx])
            {
                continue;
            }

            auto&& segment = vessels.segments[segmentIdx];
            auto p1Idx = segment.front();
            auto p2Idx = segment.back();

            for (auto pIdx : std::vector< VesselNetwork::IndexType >{{p1Idx, p2Idx}})
            {
                if (!handledPoint[pIdx])
                {
                    handledPoint[pIdx] = true;
                    auto degree = mapLargerToMax ? std::min(vessels.degrees[pIdx], maxDegree) : vessels.degrees[pIdx];
                    branchPointCountByDegree[degree]++;

                    // is it a branch point?
                    if ((degree > 2) && (degree <= maxDegree))
                    {
                        // Keep a map between the real point index and the branch point index used in the
                        // branchPointData struct
                        localBranchPointData.idxMap[pIdx] = localBranchPointData.degrees.size();

                        // Degrees
                        localBranchPointData.degrees.push_back(degree);

                        // To get the radius, we need to check all the radii of all adjacent segments
                        Real maxRadius = 0.0;
                        for (auto adjacentSegIdx : vessels.adjacency[pIdx])
                        {
                            maxRadius = std::max(maxRadius, localSegmentData.radiiSingle[adjacentSegIdx]);
                        }

                        // Diameter and capillary
                        localBranchPointData.diameter.push_back(maxRadius * 2.0);
                        bool isACapillary = isCapillary(maxRadius);
                        localBranchPointData.isCapillary.push_back(isACapillary);

                        // Count for density
                        localBranchPointData.branchPointDensity += 1.0;
                        localBranchPointData.branchPointDegree3Density += (degree == 3) ? 1.0 : 0.0;
                        localBranchPointData.branchPointDegree4Density += (degree == 4) ? 1.0 : 0.0;

                        // Count and diff between capillary and non-capillary
                        localBranchPointData.branchPointDensityCapillary += (isACapillary) ? 1.0 : 0.0;
                        localBranchPointData.branchPointDegree3DensityCapillary +=
                            (isACapillary && (degree == 3)) ? 1.0 : 0.0;
                        localBranchPointData.branchPointDegree4DensityCapillary +=
                            (isACapillary && (degree == 4)) ? 1.0 : 0.0;

                        localBranchPointData.branchPointDensityNonCapillary += (!isACapillary) ? 1.0 : 0.0;
                        localBranchPointData.branchPointDegree3DensityNonCapillary +=
                            (!isACapillary && (degree == 3)) ? 1.0 : 0.0;
                        localBranchPointData.branchPointDegree4DensityNonCapillary +=
                            (!isACapillary && (degree == 4)) ? 1.0 : 0.0;

                        // Also accumulate the degrees into the mean values
                        localBranchPointData.degreeMean += degree;
                        localBranchPointData.degreeMeanCapillary += (isACapillary) ? degree : 0.0;
                        localBranchPointData.degreeMeanNonCapillary += (!isACapillary) ? degree : 0.0;
                    }
                }
            }
        }

        // Branch Degree Means
        // NOTE: the densities are counters and not yet scaled for volume.
        localBranchPointData.degreeMean /= localBranchPointData.branchPointDensity;
        localBranchPointData.degreeMeanCapillary /= localBranchPointData.branchPointDensityCapillary;
        localBranchPointData.degreeMeanNonCapillary /= localBranchPointData.branchPointDensityNonCapillary;

        // Relative to volume:
        localBranchPointData.branchPointDensity /= networkCutDomainVolume;
        localBranchPointData.branchPointDegree3Density =
            (localBranchPointData.branchPointDegree3Density > 0.0)
                ? localBranchPointData.branchPointDegree3Density / networkCutDomainVolume
                : 0.0;
        localBranchPointData.branchPointDegree4Density =
            (localBranchPointData.branchPointDegree4Density > 0.0)
                ? localBranchPointData.branchPointDegree4Density / networkCutDomainVolume
                : 0.0;

        localBranchPointData.branchPointDensityCapillary /= networkCutDomainVolume;
        localBranchPointData.branchPointDegree3DensityCapillary =
            (localBranchPointData.branchPointDegree3DensityCapillary > 0.0)
                ? localBranchPointData.branchPointDegree3DensityCapillary / networkCutDomainVolume
                : 0.0;
        localBranchPointData.branchPointDegree4DensityCapillary =
            (localBranchPointData.branchPointDegree4DensityCapillary > 0.0)
                ? localBranchPointData.branchPointDegree4DensityCapillary / networkCutDomainVolume
                : 0.0;

        localBranchPointData.branchPointDensityNonCapillary /= networkCutDomainVolume;
        localBranchPointData.branchPointDegree3DensityNonCapillary =
            (localBranchPointData.branchPointDegree3DensityNonCapillary > 0.0)
                ? localBranchPointData.branchPointDegree3DensityNonCapillary / networkCutDomainVolume
                : 0.0;
        localBranchPointData.branchPointDegree4DensityNonCapillary =
            (localBranchPointData.branchPointDegree4DensityNonCapillary > 0.0)
                ? localBranchPointData.branchPointDegree4DensityNonCapillary / networkCutDomainVolume
                : 0.0;

        LogI << "Branch Point Info: ";
        for (auto deg : range(3, highestDegree + 1))
        {
            if (branchPointCountByDegree[deg] != 0)
            {
                LogContinue << "d" << deg << " = " << branchPointCountByDegree[deg] << ", ";
            }
        }
        LogContinue << LogEnd;
        LogI << "Branch Degree Means - all: " << localBranchPointData.degreeMean
             << " - capillary: " << localBranchPointData.degreeMeanCapillary
             << " - non-capillary: " << localBranchPointData.degreeMeanNonCapillary << LogEnd;

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 5 - Analyze volume
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (!useVVF)
        {
            LogI << "Skipping calculating vessel volume fractions." << LogEnd;
        }
        else
        {
            auto cutVolumeVVFBB = cutBoundingBox(maskVolume->boundingBox, volumeCutXYVVF, volumeCutZVVF);
            LogD << "Cut Volume Bounding Box: [ ( " << cutVolumeVVFBB.first[0] << ", " << cutVolumeVVFBB.first[1]
                 << ", " << cutVolumeVVFBB.first[2] << " ), ( " << cutVolumeVVFBB.second[0] << ", "
                 << cutVolumeVVFBB.second[1] << ", " << cutVolumeVVFBB.second[2] << " ) ]" << LogEnd;
            auto cutVolumeVVFLocalizedBB =
                cutBoundingBox(maskVolume->boundingBox, volumeCutXYVVFLocalized, volumeCutZVVFLocalized);
            LogD << "Cut Volume Bounding Box for localized evaluation: [ ( " << cutVolumeVVFLocalizedBB.first[0] << ", "
                 << cutVolumeVVFLocalizedBB.first[1] << ", " << cutVolumeVVFLocalizedBB.first[2] << " ), ( "
                 << cutVolumeVVFLocalizedBB.second[0] << ", " << cutVolumeVVFLocalizedBB.second[1] << ", "
                 << cutVolumeVVFLocalizedBB.second[2] << " ) ]" << LogEnd;

            LogI << "Calculating vessel volume fractions." << LogEnd;
            // This is mostly trivial.
            // Iterate all voxels and count using the given mask
            size_t vesselVoxelCount = 0;
            auto volumeVVFCylinderVoxelCount =
                processVoxels(cutVolumeVVFBB, isCylindric, [&](auto&& x, auto&& y, auto&& z) {
                    auto mask = maskVolume->data[maskVolume->index(x, y, z)];
                    vesselVoxelCount += (mask == 0) ? 0 : 1;
                });

            // This is mostly trivial. Do for another BB for localized VVF info
            // Iterate all voxels and count using the given mask
            size_t vesselVoxelLocalCount = 0;
            auto volumeVVFLocalCylinderVoxelCount =
                processVoxels(cutVolumeVVFLocalizedBB, isCylindric, [&](auto&& x, auto&& y, auto&& z) {
                    auto mask = maskVolume->data[maskVolume->index(x, y, z)];
                    vesselVoxelLocalCount += (mask == 0) ? 0 : 1;
                });

            // The fraction:
            globalVolumeData.vesselVolumeFraction =
                static_cast< Real >(vesselVoxelCount) / static_cast< Real >(volumeVVFCylinderVoxelCount);

            // LOCAL !

            // We know the percentage of capillaries -> so we derive the fractions for capillaries and non-capillaries
            // accordingly: Keep in mind: this is not perfectly accurate ... but suffices for the statistics.
            globalVolumeData.vesselVolumeFractionLocal =
                static_cast< Real >(vesselVoxelLocalCount) / static_cast< Real >(volumeVVFLocalCylinderVoxelCount);

            auto segmentVolume = vessels.volumeCapillaries + vessels.volumeNonCapillaries;
            auto capFrac = vessels.volumeCapillaries / segmentVolume;
            auto nonCapFrac = vessels.volumeNonCapillaries / segmentVolume;

            globalVolumeData.vesselVolumeFractionLocalCapillary = globalVolumeData.vesselVolumeFractionLocal * capFrac;
            globalVolumeData.vesselVolumeFractionLocalNonCapillary =
                globalVolumeData.vesselVolumeFractionLocal * nonCapFrac;

            // I nearly forgot it -> we want it to be percentages
            globalVolumeData.vesselVolumeFraction *= 100.0;
            globalVolumeData.vesselVolumeFractionLocal *= 100.0;
            globalVolumeData.vesselVolumeFractionLocalCapillary *= 100.0;
            globalVolumeData.vesselVolumeFractionLocalNonCapillary *= 100.0;

            LogD << "Volume contains " << vesselVoxelCount << " vessel voxels in a total of "
                 << volumeVVFCylinderVoxelCount << " voxels." << LogEnd;
            LogD << "Network-local Volume contains " << vesselVoxelLocalCount << " vessel voxels in a total of "
                 << volumeVVFLocalCylinderVoxelCount << " voxels." << LogEnd;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 6 - Analyze distance transformed volume
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if (!useEVD)
        {
            LogI << "Skipping calculating EVD fractions." << LogEnd;
        }
        else
        {
            auto cutVolumeEVDBB = cutBoundingBox(maskVolume->boundingBox, volumeCutXYEVD, volumeCutZEVD);
            LogD << "Cut Volume Bounding Box for EVD: [ ( " << cutVolumeEVDBB.first[0] << ", "
                 << cutVolumeEVDBB.first[1] << ", " << cutVolumeEVDBB.first[2] << " ), ( " << cutVolumeEVDBB.second[0]
                 << ", " << cutVolumeEVDBB.second[1] << ", " << cutVolumeEVDBB.second[2] << " ) ]" << LogEnd;

            LogI << "Calculating EVD fractions." << LogEnd;

            // This is mostly trivial.
            // Iterate all voxels and count using the given mask
            size_t evdVoxelsUsed = 0;
            size_t evdIgnoredByHistogram = 0;
            size_t evdIgnoredByMinimum = 0;
            size_t evdIgnoredByMaximum = 0;
            size_t evdIgnoredByVolumeCrop = 0;

            Real evdHistBinWidth = (evdHistMax - evdHistMin) / static_cast< Real >(evdNumBins);

            // Init the bin container
            globalEVDData.histBinsCenters.resize(evdNumBins, 0.0);
            size_t binCenterNum = 0;
            std::generate(globalEVDData.histBinsCenters.begin(), globalEVDData.histBinsCenters.end(), [&]() {
                Real center =
                    evdHistMin + (evdHistBinWidth / 2.0) + (static_cast< Real >(binCenterNum) * evdHistBinWidth);
                binCenterNum++;
                return center;
            });

            std::vector< size_t > evdHistogram(evdNumBins, 0);
            for (auto z : range(dtVolume->sizeZ))
            {
                for (auto y : range(dtVolume->sizeY))
                {
                    for (auto x : range(dtVolume->sizeX))
                    {
                        // Get the current value
                        Real evd = evdScale * dtVolume->data[dtVolume->index(x, y, z)];
                        evdVoxelsUsed++;

                        auto p = Vec3({{static_cast< Real >(x), static_cast< Real >(y), static_cast< Real >(z)}});
                        if (!pointInVolume(p, cutVolumeEVDBB, isCylindric))
                        {
                            evdIgnoredByVolumeCrop++;
                            continue;
                        }
                        if ((evd < minEVD) && useMinEVD)
                        {
                            evdIgnoredByMinimum++;
                            continue;
                        }
                        if ((evd > maxEVD) && useMaxEVD)
                        {
                            evdIgnoredByMaximum++;
                            continue;
                        }

                        globalEVDData.max = std::max(evd, globalEVDData.max);
                        globalEVDData.mean += evd;

                        // Which bin?
                        auto binNum = calcBinNum(evd, evdNumBins, evdHistMin, evdHistMax);

                        // Inside the histogram?
                        if ((binNum < 0) || (binNum >= evdNumBins))
                        {
                            evdIgnoredByHistogram++;
                            continue;
                        }

                        evdHistogram[binNum]++;
                    }
                }
            }

            // As we want relative amounts, scale by the used voxels and provide some debug output
            LogD << "Histogram - ignored of " << evdVoxelsUsed << " voxels: " << LogEnd;
            LogD << " - by volume crop: " << evdIgnoredByVolumeCrop << LogEnd;
            LogD << " - by defined min: " << evdIgnoredByMinimum << LogEnd;
            LogD << " - by defined max: " << evdIgnoredByMaximum << LogEnd;
            LogD << " - by histogram boundaries: " << evdIgnoredByHistogram << LogEnd;

            size_t evdReallyUsed = evdVoxelsUsed - evdIgnoredByHistogram - evdIgnoredByMinimum - evdIgnoredByMaximum -
                                   evdIgnoredByVolumeCrop;
            std::vector< Real > evdHistogramPerc(evdNumBins, 0.0);
            for (auto binNum : range(evdHistogram.size()))
            {
                // get the bin range:
                auto binLeft = evdHistMin + binNum * evdHistBinWidth;
                auto binRight = evdHistMin + (binNum + 1) * evdHistBinWidth;
                auto perc = 100.0 * static_cast< Real >(evdHistogram[binNum]) / static_cast< Real >(evdReallyUsed);

                LogD << "[ " << binLeft << ", " << binRight << " )"
                     << " - " << perc << " ( " << evdHistogram[binNum] << " )" << LogEnd;
                evdHistogramPerc[binNum] = perc;
            }

            globalEVDData.mean /= static_cast< Real >(evdVoxelsUsed);
            LogI << "Calculated EVD histogram. Mean is " << globalEVDData.mean << ", max is " << globalEVDData.max
                 << "." << LogEnd;

            globalEVDData.histMin = evdHistMin;
            globalEVDData.histMax = evdHistMax;
            globalEVDData.histNumBins = evdNumBins;

            globalEVDData.discretized = evdHistogram;
            globalEVDData.discretizedRelative = evdHistogramPerc;
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 7 - Writing local data to files
        //
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        LogI << "Writing data ..." << LogEnd;
        write(outputDir, localBranchPointData);
        write(outputDir, localSegmentData);
        write(outputDir, vessels, vessels.points, localSegmentData, localBranchPointData);
        if (useVVF)
        {
            write(outputDir, globalVolumeData);
        }
        if (useEVD)
        {
            write(outputDir, globalEVDData);
        }

        LogI << "Done ..." << LogEnd;
    }
} // namespace nogo
