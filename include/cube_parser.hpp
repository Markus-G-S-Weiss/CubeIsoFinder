/*
 * CubeIsoFinder
 * File: cube_parser.hpp
 *
 * Description:
 *   Contains declarations for cube file parsing, unit conversion, and integration functions.
 *
 * Author: Markus G. S. Weiss
 * Created: 2025-02-13
 *
 * License: GNU GPL v3.0
 *   See LICENSE file in the project root for full license information.
 */

#ifndef CUBE_PARSER_HPP
#define CUBE_PARSER_HPP

#include <string>
#include <vector>

// Structure representing the header information of a cube file.
// ----- Cube File Parsing Structures -----
//
// CubeHeader holds information about the cube file. It contains the
// first two comment lines, number of atoms, the origin, grid dimensions,
// axis vectors (each with the voxel count and 3 vector components),
// the calculation type (e.g., "Q-Chem", "ORCA", "Generic"),
// and a flag indicating whether the data are orbital (true) or density (false).
struct CubeHeader {
    std::string comment1;
    std::string comment2;
    int numAtoms;
    double origin[3];
    int dims[3];              // Number of voxels in x, y, and z directions.
    double axisVectors[3][4]; // Each row: [n, ax, ay, az] for the axis (n is the voxel count).
    std::string calcType;     // "Q-Chem", "ORCA", or "Generic".
    bool isOrbital;           // True if orbital data; false if density data.
};

// CubeData holds a CubeHeader and a flat vector of doubles that contains
// the volumetric grid data in the order it was read.
struct CubeData {
    CubeHeader header;
    std::vector<double> values;
};

// Helper function declarations.
std::string trim(const std::string &s);
bool icontains(const std::string &data, const std::string &substr);

// Cube file parsing functions.
CubeData readCubeFile(const std::string &filename);
double computeVoxelVolume(const CubeHeader &header);

// Unit detection and conversion functions.
bool detectAngstrom(const CubeHeader &header);
double convertDensity(double nativeDensity, bool nativeIsAngstrom);
double convertOrbital(double nativeOrbital, bool nativeIsAngstrom);

// Integration functions for density data.
double computeIsovalueFromPercentage_Density(const std::vector<double> &values, double percent, bool positive);
double computePercentageFromIsovalue_Density(const std::vector<double> &values, double isovalue, bool positive);

// Integration functions for orbital data.
double computeIsovalueFromPercentage_Orbital(const std::vector<double> &values, double percent, bool positive);
double computePercentageFromIsovalue_Orbital(const std::vector<double> &values, double isovalue, bool positive);

#endif // CUBE_PARSER_HPP

