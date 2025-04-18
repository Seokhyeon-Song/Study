cmake_minimum_required(VERSION 3.20)

project(RMPS LANGUAGES CXX)

# Enable C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MKL CONFIG REQUIRED)
message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(RMPS ${SOURCES})

target_compile_options(RMPS PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(RMPS PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(RMPS PUBLIC $<LINK_ONLY:MKL::MKL>)