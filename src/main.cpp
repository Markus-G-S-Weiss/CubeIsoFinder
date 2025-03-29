/*
 * CubeIsoFinder
 * File: main.cpp
 *
 * Description:
 *   Entry point for the CubeIsoFinder application. Handles command-line parsing
 *   and calls functions from cube_parser.hpp to process cube files.
 *
 * Author: Markus G. S. Weiss
 * Created: 2025-02-13
 *
 * License: GNU GPL v3.0
 *   See LICENSE file in the project root for full license information.
 */

#include "cube_parser.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

// Print usage information.
void printUsage(const char *progName) {
    std::cerr << "Usage:\n"
              << "  " << progName << " <cube_file> (-p <percentage> | -v <isovalue>) [-s pos|neg]\n\n"
              << "Options:\n"
              << "  -p <percentage>   Compute the isovalue corresponding to the given percentage of charge.\n"
              << "  -v <isovalue>     Compute the percentage of total charge enclosed by the given isovalue.\n"
              << "  -s pos|neg        (For density files) Choose positive (default) or negative values for integration.\n";
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printUsage(argv[0]);
        return 1;
    }

    std::string cubeFilename;
    bool usePercentage = false;
    bool useIsovalue = false;
    double inputValue = 0.0;
    bool positive = true; // Default for density data.

    cubeFilename = argv[1];

    // Process command-line arguments.
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-p" && i + 1 < argc) {
            usePercentage = true;
            inputValue = std::stod(argv[++i]);
        }
        else if (arg == "-v" && i + 1 < argc) {
            useIsovalue = true;
            inputValue = std::stod(argv[++i]);
        }
        else if (arg == "-s" && i + 1 < argc) {
            std::string signStr = argv[++i];
            if (signStr == "pos")
                positive = true;
            else if (signStr == "neg")
                positive = false;
            else {
                std::cerr << "Error: Invalid sign option. Use 'pos' or 'neg'.\n";
                return 1;
            }
        }
        else {
            std::cerr << "Error: Unknown option " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Exactly one of -p or -v must be specified.
    if (usePercentage == useIsovalue) {
        std::cerr << "Error: You must specify exactly one of -p (percentage) or -v (isovalue).\n";
        printUsage(argv[0]);
        return 1;
    }

    try {
        // Read the cube file.
        CubeData cube = readCubeFile(cubeFilename);
        // Compute voxel volume using the grid's axis vectors.
        double voxelVolume = computeVoxelVolume(cube.header);
        // Determine the native unit.
        bool nativeIsAngstrom = detectAngstrom(cube.header);
        std::string nativeUnit = nativeIsAngstrom ? "Å" : "bohr";
        std::string convUnit = nativeIsAngstrom ? "bohr" : "Å";

        std::cout << "Processing file: " << cubeFilename << "\n";
        std::cout << "Calculation type detected: " << cube.header.calcType << "\n";
        std::cout << "Data type detected: " << (cube.header.isOrbital ? "Orbital" : "Density") << "\n";
        std::cout << "Grid dimensions: " << cube.header.dims[0] << " x "
                  << cube.header.dims[1] << " x " << cube.header.dims[2] << "\n";
        std::cout << "Voxel volume: " << voxelVolume << " " << nativeUnit << "^3\n";

        // Compute the total integrated density.
        double totalIntegrated = 0.0;
        if (cube.header.isOrbital) {
            // For orbital data, integrate the square (orbital density).
            for (double v : cube.values)
                totalIntegrated += v * v;
            totalIntegrated *= voxelVolume;
            std::cout << "Total integrated orbital density: " << totalIntegrated << "\n";
        }
        else {
            // For density data, integrate the values directly.
            for (double v : cube.values)
                totalIntegrated += v;
            totalIntegrated *= voxelVolume;
            std::cout << "Total integrated electron density: " << totalIntegrated << "\n";
        }

        // For orbital data, if sign is not explicitly chosen, decide automatically.
        if (cube.header.isOrbital) {
            double posTotal = 0.0, negTotal = 0.0;
            for (double v : cube.values) {
                if (v > 0)
                    posTotal += v * v;
                else if (v < 0)
                    negTotal += v * v;
            }
            if (!((positive && posTotal > 0) || (!positive && negTotal > 0))) {
                positive = (posTotal >= std::abs(negTotal));
            }
        }

        // Depending on whether a percentage or a specific isovalue was provided, compute the mapping.
        if (usePercentage) {
            if (cube.header.isOrbital) {
                std::cout << "Integrating (in orbital mode) to reach " << inputValue << "% of the total quantity...\n";
                double isovalue_native = computeIsovalueFromPercentage_Orbital(cube.values, inputValue, positive);
                double isovalue_converted = convertOrbital(isovalue_native, nativeIsAngstrom);
                std::cout << "Isovalue (orbital) corresponding to " << inputValue << "%:\n"
                          << "  " << isovalue_native << " (native, electrons/" << nativeUnit << "^(3/2))\n"
                          << "  " << isovalue_converted << " (converted, electrons/" << convUnit << "^(3/2))\n";

                double integratedAbove = 0.0;
                for (double v : cube.values) {
                    if (v * v >= isovalue_native * isovalue_native)
                        integratedAbove += v * v;
                }
                integratedAbove *= voxelVolume;
                std::cout << "Integrated orbital density above threshold (native): " << integratedAbove << "\n";
                double enclosedPercentage = computePercentageFromIsovalue_Orbital(cube.values, isovalue_native, positive);
                std::cout << "Computed percentage of total orbital density above threshold: "
                          << enclosedPercentage << "%\n";
            }
            else {
                std::cout << "Integrating (in density mode) to reach " << inputValue << "% of the total quantity...\n";
                double isovalue_native = computeIsovalueFromPercentage_Density(cube.values, inputValue, positive);
                double isovalue_converted = convertDensity(isovalue_native, nativeIsAngstrom);
                std::cout << "Isovalue (density) corresponding to " << inputValue << "%:\n"
                          << "  " << isovalue_native << " (native, electrons/" << nativeUnit << "^3)\n"
                          << "  " << isovalue_converted << " (converted, electrons/" << convUnit << "^3)\n";

                double integratedAbove = 0.0;
                for (double v : cube.values) {
                    if ((positive && v >= isovalue_native) || (!positive && v <= isovalue_native))
                        integratedAbove += v;
                }
                integratedAbove *= voxelVolume;
                std::cout << "Integrated electron density above threshold (native): " << integratedAbove << "\n";
                double enclosedPercentage = computePercentageFromIsovalue_Density(cube.values, isovalue_native, positive);
                std::cout << "Computed percentage of total electron density above threshold: "
                          << enclosedPercentage << "%\n";
            }
        }
        else {
            if (cube.header.isOrbital) {
                double percentage = computePercentageFromIsovalue_Orbital(cube.values, inputValue, positive);
                std::cout << "For orbital data, the percentage of total charge enclosed by isovalue " 
                          << inputValue << " (electrons/" << nativeUnit << "^(3/2)) is: " 
                          << percentage << "%\n";
                std::cout << "Converted isovalue: " << convertOrbital(inputValue, nativeIsAngstrom)
                          << " electrons/" << convUnit << "^(3/2)\n";
            }
            else {
                double percentage = computePercentageFromIsovalue_Density(cube.values, inputValue, positive);
                std::cout << "For density data, the percentage of total charge enclosed by isovalue " 
                          << inputValue << " (electrons/" << nativeUnit << "^3) is: " 
                          << percentage << "%\n";
                std::cout << "Converted isovalue: " << convertDensity(inputValue, nativeIsAngstrom)
                          << " electrons/" << convUnit << "^3\n";
            }
        }
    }
    catch (const std::exception &ex) {
        std::cerr << "Exception encountered: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}

