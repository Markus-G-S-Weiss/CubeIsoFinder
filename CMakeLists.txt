cmake_minimum_required(VERSION 3.10)
project(CubeIsoFinder VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include header files from the include directory.
include_directories(${PROJECT_SOURCE_DIR}/include)

# Define the executable target.
add_executable(CubeIsoFinder src/main.cpp src/cube_parser.cpp)

