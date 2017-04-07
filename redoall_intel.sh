#!/bin/bash
# Proper header for a Bash script.

#
# Rebuild all (or cleans) in LMAT for INTEL
#
# You can provide a parameter for choosing the build profile of CMake:
# 	D for Debug (default!)
#	R for Release
#	I for RelWithDebInfo
#	M for MinSizeRel
# For just cleaning the parameter is 'clean'
#
# JMMM - April 2017 - Rel. 0.2
#

# This is specific to LC
echo "Preparing environment for INTEL compilers..."
#use ic-15.0.187 #13.1.163 #14.0.211 #15.0.187

# This is general
comp="-D COMPILER_FAMILY=intel" 
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
    cmake $comp -D CMAKE_BUILD_TYPE=Release .
  elif [ $1 = "I" ]
  then
    echo "  CMake profile will be RelWithDebInfo"
    cmake $comp -D CMAKE_BUILD_TYPE=RelWithDebInfo .
  elif [ $1 = "M" ]
  then
    echo "  CMake profile will be MinSizeRel"
    cmake $comp -D CMAKE_BUILD_TYPE=MinSizeRel .
  elif [ $1 = "D" ]
  then
    echo "  CMake profile will be Debug"
    cmake $comp -D CMAKE_BUILD_TYPE=Debug .
  fi
else
  echo "  CMake profile will be the default one"
  cmake $comp .
fi
#make -j 16
make
