cmake_minimum_required(VERSION 3.20)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project(RMPS LANGUAGES CXX)

# Enable C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(RMPS ${SOURCES})
target_include_directories(RMPS PRIVATE ${CMAKE_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/external/OpenBLAS/include)
target_link_directories(RMPS PRIVATE ${PROJECT_SOURCE_DIR}/external/OpenBLAS/lib)

target_link_libraries(RMPS PRIVATE openblas)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Using GCC compiler")
endif()