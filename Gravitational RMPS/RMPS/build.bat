@echo off
echo Cleaning previous build...
rmdir /s /q build
mkdir build
cd build

echo Running CMake...
cmake .. -G "Visual Studio 17 2022"
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