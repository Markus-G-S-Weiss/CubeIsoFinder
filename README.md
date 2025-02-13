# CubeIsoFinder

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/Markus-G-S-Weiss/CubeIsoFinder.svg)]()

*Created by Markus G. S. Weiss on 2025-02-13.*

CubeIsoFinder is a lightweight, modern C++ tool built with CMake that parses cube files, integrates density or orbital data, and computes isovalues from user-specified percentages.

## Features

- Parse cube files to extract volumetric grid data.
- Auto-detect cube file format and data type.
- Compute voxel volumes based on grid axis vectors.
- Map a specified percentage of integrated data to an isovalue (or vice versa).
- Automatic unit conversion between Angstroms and bohrs.

## Programs Included

CubeIsoFinder provides a single executable that supports two modes:
- **Percentage Mode (-p)**: Computes the isovalue corresponding to a given percentage of integrated data.
- **Isovalue Mode (-v)**: Computes the percentage of integrated data above a given isovalue.

## Installation Instructions

### Requirements

- A C++ compiler supporting C++17 or later.
- CMake (version 3.10 or later).

### Build Instructions

1. Clone the repository:

   ```
   git clone https://github.com/Markus-G-S-Weiss/CubeIsoFinder.git
   ```

2. Change directory into the project folder:

   ```
   cd CubeIsoFinder
   ```

3. Create and navigate to the build directory:

   ```
   mkdir build
   cd build
   ```

4. Configure the project with CMake:

   ```
   cmake ..
   ```

5. Build the project:

   ```
   cmake --build .
   ```

## Usage

Run the executable with the following syntax:

   ```
   ./CubeIsoFinder <cube_file> (-p <percentage> | -v <isovalue>) [-s pos|neg]
   ```

**Parameters:**

- `<cube_file>`: Path to the cube file.
- `-p <percentage>`: Compute the isovalue corresponding to the specified percentage of integrated data.
- `-v <isovalue>`: Compute the percentage of integrated data above the specified isovalue.
- `-s pos|neg`: For density data, specify positive (default) or negative integration.

## Example Usage

To compute the isovalue for 50% of the integrated data:

   ```
   ./CubeIsoFinder example.cube -p 50
   ```

To compute the percentage of integrated data above an isovalue of 0.05:

   ```
   ./CubeIsoFinder example.cube -v 0.05
   ```

## License

CubeIsoFinder is released under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## Authors

- **Markus G. S. Weiss**

  Primary developer and maintainer of CubeIsoFinder.

This project builds upon established methodologies for cube file analysis and integration, and is continually improved with contributions from the open-source community.

## Contact

For issues, suggestions, or contributions, please contact Markus G. S. Weiss or create issues or pull requests.

