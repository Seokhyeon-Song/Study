@echo off
setlocal

set VCPKG_ROOT=C:/Users/hypercube256/Development/vcpkg
set VCPKG_TRIPLET=x64-mingw-dynamic
set VCPKG_LIBDIR=%VCPKG_ROOT%/installed/%VCPKG_TRIPLET%/lib

echo Cleaning previous build...
rmdir /s /q build 2>nul
mkdir build
cd build

echo Running CMake...
cmake .. -G "MinGW Makefiles" ^
    -DCMAKE_BUILD_TYPE=Release

if %errorlevel% neq 0 (
    echo CMake configuration failed!
    exit /b %errorlevel%
)

echo Building project...
cmake --build . --config Release
if %errorlevel% neq 0 (
    echo Build failed!
    exit /b %errorlevel%
)

echo Build complete!
endlocal