#!/bin/bash
# Proper header for a Bash script.

#
# Rebuild all (or cleans) in LMAT for CLANG (with OpenMP support!)
#
# You can provide a parameter for choosing the build profile of CMake:
# 	D for Debug (default!)
#	R for Release
#	I for RelWithDebInfo
#	M for MinSizeRel
# For just cleaning the parameter is 'clean'
#
# JMMM - April 2017 - Rel. 0.3

echo "Preparing environment for CLANG..."
cmakebin="cmake"
comp="-D COMPILER_FAMILY=clang" 
echo "Performing the total reconf and rebuild..."
make clean
rm -Rf *.cmake Makefile src/Makefile src/CMakeCache.txt src/*.cmake CMakeCache.txt src/kmerdb/all_headers.hpp CMakeFiles/ src/CMakeFiles/
find bin/ -type f  ! -name "*.*"  -delete
find third-party/. -mindepth 1 ! -name "README" -exec rm -rf "{}" +
if [ -n "$1" ]
then
  if [ $1 = "clean" ]
  then
    echo "  Just cleaning, nothing built!"
    exit 0
  elif [ $1 = "R" ]
  then
    echo "  CMake profile will be Release"
    $cmakebin $comp -D CMAKE_BUILD_TYPE=Release .
  elif [ $1 = "I" ]
  then
    echo "  CMake profile will be RelWithDebInfo"
    $cmakebin $comp -D CMAKE_BUILD_TYPE=RelWithDebInfo .
  elif [ $1 = "M" ]
  then
    echo "  CMake profile will be MinSizeRel"
    $cmakebin $comp -D CMAKE_BUILD_TYPE=MinSizeRel .
  elif [ $1 = "D" ]
  then
    echo "  CMake profile will be Debug"
    $cmakebin $comp -D CMAKE_BUILD_TYPE=Debug .
  fi
else
  echo "  CMake profile will be the default one"
  $cmakebin $comp .
fi
make 
