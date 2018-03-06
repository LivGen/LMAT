#!/bin/bash

# redoall - JMMM - March 2018 - Rel. 0.1

defaultcomp="gnu" # Default compiler

print_usage(){
echo " "
echo "redoall: Easily rebuild all (or cleans) LMAT with CMake"
echo " "
echo "usage: redoall [profile] [compiler]"
echo " "
echo "The 1st optional parameter chooses the build profile of CMake:"
echo "  D for Debug"
echo "  R for Release (default!)"
echo "  I for RelWithDebInfo"
echo "  M for MinSizeRel"
echo "For just cleaning the parameter is 'clean'"
echo " "
echo "The 2nd optional parameter selects the compiler family:"
echo "  gnu for using GCC"
echo "  intel for using Intel compilers"
echo "  clang for using clang compilers"
echo "  ibmpwr9 for compiling in Power 9 with IBM compilers"
echo "  ibmpwr8 for compiling in Power 8 with IBM compilers"
echo "Default compiler : $defaultcomp" 
exit 1
}

if [ "$#" -eq 2 ] && ([ $2 = "gnu" ] || [ $2 = "intel" ] || [ $2 = "clang" ] || [ $2 = "ibmpwr9" ] || [ $2 = "ibmpwr8" ]); then
    echo "Preparing environment for $2 compilers..."
    comp="-D COMPILER_FAMILY=$2"
elif [ "$#" -eq 2 ]; then
    echo "ERROR! Unsupported compiler"
    print_usage
elif [ "$#" -lt 2 ]; then
    echo "Using default compiler: $defaultcomp"
    comp="-D COMPILER_FAMILY=$defaultcomp"
else
    echo "ERROR! Invalid number of parameters"
    print_usage
fi
 
if hash cmake3 2>/dev/null; then
    cmakebin="cmake3"
elif hash cmake 2>/dev/null; then
    cmakebin="cmake"
else
    echo "ERROR! CMake not found in the system!"
    echo "NOTE: CMake >= 3.2 required to build LMAT"
fi
echo "Using to build LMAT the following CMake release:"
$cmakebin --version | head -1

echo "Performing the total reconf and rebuild..."
make clean 2>/dev/null
rm -Rf *.cmake Makefile src/Makefile src/CMakeCache.txt src/*.cmake CMakeCache.txt src/kmerdb/all_headers.hpp CMakeFiles/ src/CMakeFiles/
find bin/ -type f  ! -name "*.*"  -delete
find third-party/. -mindepth 1 ! -name "README" -exec rm -rf "{}" +

if [ -n "$1" ]; then
    if [ $1 = "clean" ]; then
        echo "  Just cleaning, nothing built!"
        exit 0
    elif [ $1 = "R" ]; then
        echo "  CMake profile will be Release"
        $cmakebin $comp -D CMAKE_BUILD_TYPE=Release .
    elif [ $1 = "I" ]; then
        echo "  CMake profile will be RelWithDebInfo"
        $cmakebin $comp -D CMAKE_BUILD_TYPE=RelWithDebInfo .
    elif [ $1 = "M" ]; then
        echo "  CMake profile will be MinSizeRel"
        $cmakebin $comp -D CMAKE_BUILD_TYPE=MinSizeRel .
    elif [ $1 = "D" ]; then
        echo "  CMake profile will be Debug"
        $cmakebin $comp -D CMAKE_BUILD_TYPE=Debug .
    else
	echo "ERROR! Unsupported build type"
        print_usage
    fi
else
    echo "Using default value for build"
    $cmakebin $comp .
fi
make 
