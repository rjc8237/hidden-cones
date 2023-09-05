cmake_minimum_required(VERSION 3.20)
project(ParametrizeMesh VERSION 0.0.1)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# Get required packages
list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(eigen)
include(spdlog)
include(libigl)

# Install executables to bin directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Build cone generation library
add_library(ConeGenerationLib
  src/cone_energy.cpp
  src/cone_validity.cpp
  src/generate_flat_cones.cpp
  src/generate_optimal_cones.cpp
)
target_include_directories(ConeGenerationLib PUBLIC src)
target_link_libraries(ConeGenerationLib PUBLIC
  igl::core
  Eigen3::Eigen
  spdlog::spdlog
)

