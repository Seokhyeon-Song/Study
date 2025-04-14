cmake_minimum_required(VERSION 3.20)

project(RMPS LANGUAGES CXX)

# Enable C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

find_package(MKL REQUIRED)

add_executable(RMPS ${SOURCES})

target_link_libraries(RMPS MKL::MKL)