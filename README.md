# Vessel Network Anayzer

A software tool to extract statistically relevant features from vessel network files and their respective binary mask volumes. It is a
part of the processing pipeline used in the [Nature Protocols](https://www.nature.com/nprot) publication *'Hierarchical imaging and computational analysis of three-dimensional vascular network architecture in the entire postnatal and adult mouse brain'*.

The source has been extracted from a vast software package and was made publicly available under GPL v3. Required algorithms, tools and data structures have been header-ified for simple build and integration.

## Purpose

WIP

## Requirements

### Hardware

* A decent machine with at least 8GB RAM.
* **Tested on** an Intel i7-10875H, 16GB RAM running Debian Sid Linux.

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
## Usage

WIP

## Example

WIP

## License

This software is licensed under [GNU General Public License v3.0.](LICENSE.txt)

This software utilizes other components under different licenses. Their respective license information can be found at:

* https://github.com/ukoethe/vigra
* https://github.com/kazuho/picojson

## Authors

* Sebastian Eichelbaum: [www.nemtics.com](http://www.nemtics.com)


