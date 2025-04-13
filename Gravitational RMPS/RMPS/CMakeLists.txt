cmake_minimum_required(VERSION 3.20)
project(RMPS)

# Enable C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen via vcpkg (no need to install separately if using vcpkg)
set(CMAKE_TOOLCHAIN_FILE "C:/Users/hypercube256/Development/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "vcpkg toolchain file")
set(CMAKE_PREFIX_PATH "C:/Users/hypercube256/Development/vcpkg/installed/x64-mingw-dynamic")
find_package(Eigen3 REQUIRED)

# Add your source files here
add_executable(RMPS src/main.cpp)
target_include_directories(RMPS PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link Eigen (header-only) and system math libraries for BLAS/LAPACK
target_link_libraries(RMPS PRIVATE Eigen3::Eigen)