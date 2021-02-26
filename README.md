# Vessel Network Anayzer

A software tool to extract statistically relevant features from vessel network
files and their respective binary mask volumes.  It is a part of the processing
pipeline used in the [Nature Protocols](https://www.nature.com/nprot)
publication *'Hierarchical imaging and computational analysis of
three-dimensional vascular network architecture in the entire postnatal and
adult mouse brain'*.

It collects information on bifurcations, segments, vascular volume and
extra-vascular distance. The information is stored as CSV files per subject and
can be used to create the plots as shown in the paper. This tool does *not*
create the plots nor the 3D renderings shown in the paper. The plots have been
created using [Matplotlib](https://matplotlib.org) and the 3D renderings where
done using [OpenWalnut](http://www.openwalnut.org).

The source has been extracted from a vast software package and was made
publicly available under GPL v3. Required algorithms, tools and data structures
have been header-ified for simple build and integration.

## Requirements

### Hardware

* A decent machine with at least 8GB RAM.
* **Tested using** an Intel i7-10875H, 16GB RAM running Debian Sid Linux.

### Building Requirements

* [CMake >= 3.0](https://cmake.org/)
* A C++14 capable compiler like [gcc >= 5](https://gcc.gnu.org) or [clang > 3.4](https://clang.llvm.org)
* The required libraries are provided in source via git submodules. No further libraries are required.
* **Tested using** gcc 10.2.1 and cmake 3.18.4

## Building

These are the building instructions on Linux. Mac OS X and Windows should work in a similar fashion but have not been tested.

1. Clone the project and get all the required submodules:
    ```sh
    git clone --recurse-submodules https://github.com/sebastian-eichelbaum/VesselNetworkAnayzer.git
    ```
1. Build
    ```sh
    cd VesselNetworkAnayzer
    mkdir -p build && cd build
    cmake ..
    make
    ```

## Usage and Configuration

```sh
    analyzer OUTPUT_PATH NETWORK_FILE [VOLUME_FILE] [VOLUMEDT_FILE] [SETTINGS_FILE]
```

* **OUTPUT_PATH**: the relative or absolute path where to write results. Needs
  to exist.
* **NETWORK_FILE**: the relative or absolute path to the network file.
    * The network file is a binary VTK file. ASCII VTK is not supported. You
      can convert them easily using the Python VTK bindings.
    * The network is a set of VTK POINTS and LINES, representing the vessel
      network. The VTK file needs to contain either per-point radii or per-line
      radii using VTK POINT_DATA or CELL_DATA.
    * Refer to the provided example for details.
* **VOLUME_FILE**: (optional) - the relative or absolute path to the volume
  file. To skip, specify an empty string via "".
    * The mask volume designates each voxel to be either a vessel or not.
    * The volume file is a binary mask volume in RAW format. Each voxel is a byte. A value >0 is interpreted as binary true, hence being a part of a vessel.
    * The exact size needs to be specified in the SETTINGS_FILE.
    * This file can be omitted. No vascular volume and extra-vascular distance information will be calculated in that case.
* **VOLUMEDT_FILE**: (optional) - the relative or absolute path to the distance
  transformed volume file used as cache. To skip, specify an empty string via "".
    * The software calculates the distance transformed volume of the given mask volume and can store it in this file.
    * It will be used next time to speed up calculations.
    * If omitted, the distance transform will be re-done each run.
* **SETTINGS_FILE**: (optional) - the configuration file to load. To skip, specify an
  empty string via ""
    * Allows to configure the analyzer and fine-tune some parameters in JSON format.
    * If omitted, the defaults, as shown in [config_defaults.json](config_defaults.json) will be used.
    * You should not omit it if your volume dimensions differ from the defaults.

### Configuration

This shows the default settings. This is a commented version of
[config_defaults.json](config_defaults.json).

The settings we have used for the ROI-based
evaluation of subjects and the full-brain data can be found in the example
folder for reference and differ from these defaults.

```json
{
    "//": "The volume size. The software only supports cubes.",
    "volumeSize" : 1024,

    "//": "The size of a voxel in micormeters. This value is required",
    "//": "for proper EVD calculations.",
    "pixelSize" : 0.73,

    "//": "If true, the network data is in voxel space and will be scaled",
    "//": "by pixel size. This affects segment lengths. Radii are used as is.",
    "networkIsVoxelSpace" : false,

    "//": "If true, the volume will be seen as a cylinder instead of a",
    "//": "cube. The software crops the volume to match the cylindric shape",
    "//": "before doing any volume calculations.",
    "cylindric" : false,

    "//": "If true, the whole network is assumed to be cropped already.",
    "//": "It will be used as is but the volume will be assumed to be defined",
    "//": "by networkCutXY and networkCutZ (i.e. for densities).",
    "preCroppedNetwork" : true,

    "//": "If pre-cropped is false, these values define the percentage of",
    "//": "the volume used. Cropping will be done equidistantly from the sides.",
    "networkCutXY" : 50.0,
    "networkCutZ" : 75.0,
    "//": "Like the network data, the masked volume can be cut in a similar way.",
    "//": "This defines the amount of volume to use for calculations, cropped",
    "//": "equidistantly from the sides.",
    "volumeCutXYVVF" : 95.0,
    "volumeCutZVVF" : 100.0,
    "//": "The same for the EVD calculations.",
    "volumeCutXYEVD" : 100.0,
    "volumeCutZEVD" : 100.0,

    "//": "Restrict the max degree of branch points. Useful to filter some",
    "//": "noisy branch points.",
    "restrictDegree" : true,
    "//": "If restricted, map larger branch point degree to the max?",
    "mapLargerToMax" : false,
    "//": "If restricted, use this a s max branch point degree.",
    "maxBranchpointDegree" : 6,

    "//": "Should very close branch points be merged into one?",
    "mergeCloseBranchPoints" : true,
    "//": "If merging close branch points, use this as a max distance to",
    "//": "merge in the same coordinate system as the network data.",
    "maxMergeDistance" : 2.0,

    "//": "We store the diameter per segment length in a discrete length-span",
    "//": "histogram. These values define the histogram.",
    "diameterPerLengthMin" : 1,
    "diameterPerLengthMax" : 70,
    "diameterPerLengthNumBins" : 14,

    "//": "Force a min and max range for EVD values.",
    "minEVD" : 0.10,
    "useMinEVD" : true,
    "maxEVD" : 200.0,
    "useMaxEVD" : true,

    "//": "Define a discretization scheme for EVD. This defines the histogram",
    "//": "used to categorize EVD values.",
    "evdHistMin" : 0,
    "evdHistMax" : 200.0,
    "evdNumBins" : 20,

    "//": "Scale the volume for evd calculation and saving. Set to 0 (default if",
    "//": "unset) to disable any scaling. When changing this value, previously",
    "//": "saved EVD volumes are probably not useable anymore.",
    "evdVolumeSize" : 0,
    "//": "When scaling the EVD volume, use this threshold to define the",
    "//": "vessel mask threshold in the scaled volume.",
    "evdScaleThreshold" : 192,
    "//": "When scaling the EVD volume, use n-order B-Spline interpolation.",
    "evdScaleSplineOrder" : 3,
    "//": "If true, the calculated EVD won't be saved",
    "dontSaveEVD" : false
}
```

### Results

* Branchpoint degree:
    * ```BranchPoints.degrees.capillary.csv, BranchPoints.degrees.csv, BranchPoints.degrees_mean.capillary.csv BranchPoints.degrees_mean.csv, BranchPoints.degrees_mean.non-capillary.csv, BranchPoints.degrees.non-capillary.csv```
    * The amount of in- and outgoing vessels per branchpoint
* Branchpoint density:
    * ```BranchPoints.density.capillary.csv, BranchPoints.density.csv, BranchPoints.density.degree3.capillary.csv, BranchPoints.density.degree3.csv, BranchPoints.density.degree3.non-capillary.csv, BranchPoints.density.degree4.capillary.csv, BranchPoints.density.degree4.csv, BranchPoints.density.degree4.non-capillary.csv, BranchPoints.density.non-capillary.csv```
    * The amount of branchpoints per volume. The volume is defined by the networkCutXY and Z settings.
* Branchpoint diameter:
    * ```BranchPoints.diameter.capillary.csv, BranchPoints.diameter.csv, BranchPoints.diameter.non-capillary.csv```
* Segments angle adjacency:
    * ```Segments.angleAdjacency.capillary.csv, Segments.angleAdjacency.csv, Segments.angleAdjacency.non-capillary.csv
    * The angles between adjacent segments.
* Segment count:
    * ```Segments.count.capillary.csv, Segments.count.csv, Segments.count.non-capillary.csv```
* Segment density:
    * ```Segments.density.capillary.csv, Segments.density.csv, Segments.density.non-capillary.csv```
    * The density is relative to the volume defined by networkCutXY and Z in the settings.
* Segment diameter:
    * ```Segments.diameter.capillary.csv, Segments.diameter.csv, Segments.diameter.non-capillary.csv, Segments.diameter_per_length_hist.csv```
    * The diameter_per_length_hist is a histogram showing the relation of length to diameter.
* Segment length:
    * ```Segments.lengths.capillary.csv, Segments.lengths.csv, Segments.lengths.non-capillary.csv```
* Segment radius:
    * ```Segments.radius.capillary.csv, Segments.radius.csv, Segments.radius.non-capillary.csv```
* Segment tortuosity:
    * ```Segments.tortuosity.capillary.csv, Segments.tortuosity.csv, Segments.tortuosity.non-capillary.csv```
    * Tortuosity is the relation between the sum of the lenths of the lines of each segment and the direct distance between the branchpoints at the start and end of the segment.
* Segment volume:
    * ```Segments.volume.capillary.csv, Segments.volume.csv, Segments.volume.non-capillary.csv```
* Extravascular Distance:
    * ```Volume.evd_hist.csv, Volume.evd_mean.csv, Volume.evd_relative_hist.csv```
    * The evd_relative_hist is the histogram with normalized values (percentages).
* Vascular volume fraction:
    * ```Volume.vvf.csv, Volume.vvf_local.capillary.csv, Volume.vvf_local.csv, Volume.vvf_local.non-capillary.csv```
    * The local VVF is the VVF of the cropped volume using ```volumeCutXYVVF``` and ```volumeCutZVVF```
* A RAW file. Thats a file with the name you specified on the command-line for ```VOLUMEDT_FILE```:
    * The distance transformed volume. The size is the same as the mask volume if the setting ```evdVolumeSize``` is 0.
* VesselSegments.stld:
    * The segments written as STL file, extended with radius, capillary classification, branchpoint classification, and branchpoint degree.
    * It is an ASCII file and contains a header describing the added attributes.

## Example

The examples folder contains one of the full-brain datasets that were part of
the calculations done in the paper. Due to the size of the volume data, it is
compressed.

```sh
cd example
mkdir -p results
tar xzf data.tar.gz
../build/analyzer results \
    data/isocortex_55.vtk \
    data/isocortex_55_1150x1150x1150.raw \
    dtCache.raw \
    config_fullbrain2.json
```

This calculates all the values and keeps a distance transformed volume as cache
for later runs. This is the way all the subjects have been processed. The
resulting CSV were loaded and pooled for plotting.

## Knows Issues

* The Existence of the result directory is never checked and the program runs
  without complains if the directory does not exist. No files are written in
  that case.

## License

This software is licensed under [GNU General Public License v3.0.](LICENSE.txt)

This software utilizes other components under different licenses. Their respective license information can be found at:

* https://github.com/ukoethe/vigra
* https://github.com/kazuho/picojson

## Authors

* Sebastian Eichelbaum: [www.nemtics.com](http://www.nemtics.com)


