/*
 * CubeIsoFinder
 * File: cube_parser.cpp
 *
 * Description:
 *   Implements functions for parsing cube files, unit conversion, and computing integration values.
 *
 * Author: Markus G. S. Weiss
 * Created: 2025-02-13
 *
 * License: GNU GPL v3.0
 *   See LICENSE file in the project root for full license information.
 */

#include "cube_parser.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ----- Helper Functions -----

// Trim whitespace from both ends of a string.
std::string trim(const std::string &s) {
    const char *ws = " \t\n\r\f\v";
    size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos)
        return "";
    size_t end = s.find_last_not_of(ws);
    return s.substr(start, end - start + 1);
}

// Case-insensitive search for a substring.
bool icontains(const std::string &data, const std::string &substr) {
    auto it = std::search(
        data.begin(), data.end(), substr.begin(), substr.end(),
        [](char ch1, char ch2) { return std::tolower(ch1) == std::tolower(ch2); });
    return (it != data.end());
}


// Read the cube file and populate a CubeData structure.
// Throws a runtime_error if the file cannot be opened or if data reading fails.
CubeData readCubeFile(const std::string &filename) {
    std::ifstream infile(filename);
    if (!infile)
        throw std::runtime_error("Error opening file: " + filename);

    CubeData cube;
    std::string line;

    // Read the first two comment lines.
    std::getline(infile, line);
    cube.header.comment1 = trim(line);
    std::getline(infile, line);
    cube.header.comment2 = trim(line);

    // Detect the calculation type based on keywords in the comment lines.
    if (icontains(cube.header.comment1, "ORCA") || icontains(cube.header.comment2, "ORCA"))
        cube.header.calcType = "ORCA";
    else if (icontains(cube.header.comment1, "Q-Chem") || icontains(cube.header.comment2, "Q-Chem"))
        cube.header.calcType = "Q-Chem";
    else
        cube.header.calcType = "Generic";

    // Detect whether the cube file contains orbital data or density data.
    if (icontains(cube.header.comment1, "MO") || icontains(cube.header.comment2, "MO") ||
        icontains(cube.header.comment1, "Orbital") || icontains(cube.header.comment2, "Orbital"))
        cube.header.isOrbital = true;
    else if (icontains(cube.header.comment1, "density") || icontains(cube.header.comment2, "density"))
        cube.header.isOrbital = false;
    else
        cube.header.isOrbital = true; // Default to orbital.

    // Read the line containing the number of atoms and the grid origin.
    std::getline(infile, line);
    std::istringstream iss(line);
    if (!(iss >> cube.header.numAtoms >> cube.header.origin[0] >> cube.header.origin[1] >> cube.header.origin[2]))
        throw std::runtime_error("Error reading number of atoms and origin.");

    // Read the three axis vectors.
    // Each of the next three lines contains the voxel count and the 3 vector components.
    for (int i = 0; i < 3; ++i) {
        std::getline(infile, line);
        std::istringstream iss_axis(line);
        if (!(iss_axis >> cube.header.dims[i] 
              >> cube.header.axisVectors[i][1] >> cube.header.axisVectors[i][2] >> cube.header.axisVectors[i][3])) {
            throw std::runtime_error("Error reading axis vector " + std::to_string(i));
        }
        // Also store the voxel count as the first element of each axis vector.
        cube.header.axisVectors[i][0] = cube.header.dims[i];
    }

    // Read the atom coordinate lines (one per atom).
    int numAtoms = std::abs(cube.header.numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        std::getline(infile, line);
        // Atom coordinates are read but not used in the integration.
    }

    // If the cube file is from an ORCA calculation, skip one extra header line (e.g., containing MO coefficients).
    if (cube.header.calcType == "ORCA") {
        std::getline(infile, line); // Skip extra line.
    }

    // Read the volumetric data.
    // The total number of grid points should equal dims[0] * dims[1] * dims[2].
    size_t totalPoints = static_cast<size_t>(cube.header.dims[0]) *
                         static_cast<size_t>(cube.header.dims[1]) *
                         static_cast<size_t>(cube.header.dims[2]);
    cube.values.reserve(totalPoints);
    double val;
    while (infile >> val) {
        cube.values.push_back(val);
    }
    if (cube.values.size() != totalPoints) {
        throw std::runtime_error("Error: Number of grid points read (" + std::to_string(cube.values.size()) +
                                 ") does not match expected (" + std::to_string(totalPoints) + ").");
    }
    return cube;
}

// Compute the voxel volume from the three axis vectors using the scalar triple product.
// The voxel volume is given by |a · (b × c)|, where a, b, c are the step vectors.
double computeVoxelVolume(const CubeHeader &header) {
    const double* a = header.axisVectors[0] + 1; // Skip the count element.
    const double* b = header.axisVectors[1] + 1;
    const double* c = header.axisVectors[2] + 1;
    double cross[3] = {
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0]
    };
    double vol = std::abs(a[0] * cross[0] + a[1] * cross[1] + a[2] * cross[2]);
    return vol;
}

// ----- Unit Detection & Conversion -----
//
// detectAngstrom attempts to determine whether the cube file’s coordinates
// are in Angstroms or in Bohr. It first searches for the keywords "angstrom" or "bohr"
// in the first two comment lines. If not found, it uses a heuristic based on the average
// length of the three axis vectors (if the average length > 2.0, assume Angstrom).
bool detectAngstrom(const CubeHeader &header) {
    if (icontains(header.comment1, "angstrom") || icontains(header.comment2, "angstrom"))
        return true;
    if (icontains(header.comment1, "bohr") || icontains(header.comment2, "bohr"))
        return false;
    double totalLength = 0.0;
    for (int i = 0; i < 3; ++i) {
        double ax = header.axisVectors[i][1],
               ay = header.axisVectors[i][2],
               az = header.axisVectors[i][3];
        double len = std::sqrt(ax * ax + ay * ay + az * az);
        totalLength += len;
    }
    double avgLength = totalLength / 3.0;
    return (avgLength > 2.0);
}

// ----- Unit Conversion Functions -----
//
// The conversion functions convert the computed isovalue from Bohr-based units to Angstrom-based units.
// For density, native values are in electrons/bohr^3 and must be converted to electrons/Å^3.
// For orbital isovalues, which are the square root of the density (units electrons/bohr^(3/2)),
// conversion uses a factor with exponent 1.5.
double convertDensity(double nativeDensity, bool nativeIsAngstrom) {
    // 1 bohr = 0.529177210544 Å, so 1 bohr^3 = (0.529177210544)^3 Å^3.
    const double bohr3_in_angstrom3 = std::pow(0.529177210544, 3.0);
    // Assume that if nativeIsAngstrom is false, then the native density is in electrons/bohr^3.
    // To convert from electrons/bohr^3 to electrons/Å^3, we divide by bohr3_in_angstrom3.
    return nativeDensity / bohr3_in_angstrom3;
}

double convertOrbital(double nativeOrbital, bool nativeIsAngstrom) {
    // For orbital isovalues, nativeOrbital is in electrons/bohr^(3/2)
    // To convert to electrons/Å^(3/2), divide by (0.529177210544)^(1.5).
    if (!nativeIsAngstrom) {
        return nativeOrbital / std::pow(0.529177210544, 1.5);
    } else {
        return nativeOrbital;
    }
}

// ----- Integration Functions -----
//
// For density data, integration is performed on the raw grid values.
// For orbital data, the integration is performed on the squared grid values (orbital density),
// and then the square root is taken at the threshold.
// These functions map a given percentage of the total integrated quantity to a threshold value (isovalue)
// and also compute the percentage from a given isovalue.

double computeIsovalueFromPercentage_Density(const std::vector<double> &values, double percent, bool positive) {
    // Filter the values by sign.
    std::vector<double> filtered;
    for (double v : values) {
        if (positive && v > 0) filtered.push_back(v);
        else if (!positive && v < 0) filtered.push_back(v);
    }
    if (filtered.empty())
        throw std::runtime_error("No grid points with the requested sign.");
    // Compute total integrated value.
    double total = 0.0;
    for (double v : filtered)
        total += v;
    double target = (percent / 100.0) * total;
    // Sort the filtered values.
    if (positive) {
        std::sort(filtered.begin(), filtered.end(), std::greater<double>());
        double integ = 0.0;
        for (double v : filtered) {
            integ += v;
            if (integ >= target)
                return v;
        }
    } else {
        std::sort(filtered.begin(), filtered.end());
        double integ = 0.0;
        for (double v : filtered) {
            integ += v;
            if (integ <= target)
                return v;
        }
    }
    return filtered.back();
}

double computePercentageFromIsovalue_Density(const std::vector<double> &values, double isovalue, bool positive) {
    double total = 0.0, integ = 0.0;
    for (double v : values) {
        if (positive && v > 0) {
            total += v;
            if (v >= isovalue)
                integ += v;
        } else if (!positive && v < 0) {
            total += v;
            if (v <= isovalue)
                integ += v;
        }
    }
    if (total == 0.0)
        throw std::runtime_error("Total charge for the requested sign is zero.");
    return (integ / total) * 100.0;
}

// For orbital data, we first square each grid value (to obtain orbital density),
// then sort and accumulate these squared values. When the cumulative sum reaches the target,
// we return the square root of the grid value (restoring the original orbital amplitude).
struct OrbitalPoint {
    double density; // v^2
    double value;   // original grid value v
    size_t index;   // original grid index in the cube file
};

double computeIsovalueFromPercentage_Orbital(const std::vector<double> &values, double percent, bool /*positive*/) {
    std::vector<OrbitalPoint> points;
    // Add all grid points for orbital data regardless of sign.
    for (size_t i = 0; i < values.size(); i++) {
        double v = values[i];
        points.push_back({v * v, v, i});
    }
    if (points.empty())
        throw std::runtime_error("No orbital grid points available.");

    double total = 0.0;
    for (const auto &p : points)
        total += p.density;
    double target = (percent / 100.0) * total;

    // Sort the points in descending order by density (i.e. squared value)
    std::sort(points.begin(), points.end(), [](const OrbitalPoint &a, const OrbitalPoint &b) {
        return a.density > b.density;
    });

    // Print sorted output in a file matching OpenCubMan's format
    /*std::ofstream debugFile("sorted_debug_finder.txt");
    if (!debugFile) {
        std::cerr << "Error: Unable to open file for writing sorted data." << std::endl;
    } else {
        debugFile << "Sorted list of grid point indices (i, |w|, d):" << std::endl;
        for (size_t i = 0; i < points.size(); ++i) {
            debugFile << i << ": index " << points[i].index
                      << " |w| = " << fabs(points[i].value)
                      << " d = " << points[i].density << std::endl;
        }
        debugFile.close();
    }*/

    // Accumulate the squared values until reaching the target fraction.
    double integ = 0.0;
    for (const auto &p : points) {
        integ += p.density;
        if (integ >= target)
            //return std::sqrt(p.density); // Return the square-rooted value.
            return p.value; // Return the square-rooted value.
    }
    //return std::sqrt(points.back().density);
    return points.back().value;
}


double computePercentageFromIsovalue_Orbital(const std::vector<double> &values, double isovalue, bool positive) {
    double thresholdDensity = isovalue * isovalue;
    double total = 0.0, integ = 0.0;
    for (double v : values) {
        double d = v * v;
        total += d;
        if (d >= thresholdDensity)
            integ += d;
        }
    if (total == 0.0)
        throw std::runtime_error("Total orbital density for the requested sign is zero.");
    return (integ / total) * 100.0;
}

