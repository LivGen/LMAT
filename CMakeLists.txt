# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

#
# Compiler family selection
#

# Choose compiler family by env var COMPILER_FAMILY
# Force compilers (do it before PROJECT!). Default are GNU.
if(COMPILER_FAMILY STREQUAL "intel")
  set(IS_INTEL TRUE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by INTEL Compilers...")
  INCLUDE(CMakeForceCompiler)
  MESSAGE("Enforcing Intel compilers!")
  set(CMAKE_C_COMPILER icc)
  set(CMAKE_CXX_COMPILER icpc)
  set(CMAKE_Fortran_COMPILER ifort)
elseif(COMPILER_FAMILY STREQUAL "ibmpwr9")
  set(IS_IBMPWR9 TRUE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by IBM Compilers for Power9...")
  INCLUDE(CMakeForceCompiler)
  MESSAGE("Enforcing IBM compilers!")
  set(CMAKE_C_COMPILER xlc)
  set(CMAKE_CXX_COMPILER xlc++)
  set(CMAKE_Fortran_COMPILER xlf)
elseif(COMPILER_FAMILY STREQUAL "clang")
  set(IS_CLANG TRUE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by CLANG w/OpenMP Compilers...")
  if(APPLE)  
    set(CMAKE_C_COMPILER clang-mp-3.8)
    set(CMAKE_CXX_COMPILER clang++-mp-3.8)
  else()
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)
  endif(APPLE)
else()
  set(IS_GNU TRUE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by GNU Compilers...")
  if(APPLE)  
    set(CMAKE_C_COMPILER gcc-mp-4.9)
    set(CMAKE_CXX_COMPILER g++-mp-4.9)
    set(CMAKE_Fortran_COMPILER gfortran-mp-4.9)  
  else()
    set(CMAKE_C_COMPILER gcc)
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_Fortran_COMPILER gfortran)  
  endif(APPLE)
endif()


#
# Parse user options and flags
#

# If the user specifies -D CMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug as default (!) prior to calling PROJECT()
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

# If the user specifies -D SORTED_DB on the command line, take its definition
# and dump it in the cache along with proper documentation, otherwise set it to True
IF(DEFINED USE_SORTED_DB)
  SET(USE_SORTED_DB ${SORTED_DB} CACHE BOOLEAN "Boolean for sort DB")
ELSE()
  SET(USE_SORTED_DB 1 CACHE BOOLEAN "Boolean for sort DB")
ENDIF()

# If the user specifies -D SDBIDX on the command line, take its definition
# and dump it in the cache along with proper documentation, otherwise set it to default
IF(DEFINED SDBIDX)
  SET(SDBIDX ${SDBIDX} CACHE STRING "Value for SDBIDX")
ELSE()
  SET(SDBIDX 2027 CACHE STRING "Value for SDBIDX")
ENDIF()

# If the user specifies -D TAXID_SIZE32 on the command line, take its definition
# and dump it in the cache along with proper documentation, otherwise set it to False
IF(DEFINED TAXID_SIZE32)
  SET(TAXID_SIZE32 ${TAXID_SIZE32} CACHE BOOLEAN "Boolean for taxid size: true=32, false=16 (default)")
ELSE()
  SET(TAXID_SIZE32 OFF CACHE BOOLEAN "Boolean for taxid size: true=32, false=16 (default)")
ENDIF()

# Taxid size related variables
IF(TAXID_SIZE32)
  SET(TAXID_SIZE 32)
  SET(TID_T uint32_t)
ELSE()
  SET(TAXID_SIZE 16)
  SET(TID_T uint16_t)
ENDIF()

# If the user specifies -D BOOST on the command line, take its definition
# and dump it in the cache along with proper documentation, otherwise set it to False
IF(DEFINED BOOST)
  SET(USE_BOOST ${BOOST} CACHE BOOLEAN "Boolean for using Boost library")
ELSE()
  SET(USE_BOOST 0 CACHE BOOLEAN "Boolean for using Boost library")
ENDIF()

# If the user specifies -D WITH_PJMALLOC on the command line, take its definition
# and dump it in the cache along with proper documentation, otherwise set it to True
IF(DEFINED WITH_PJMALLOC)
  SET(PJMALLOC ${WITH_PJMALLOC} CACHE BOOLEAN "Boolean for using PJmalloc library")
ELSE()
  SET(PJMALLOC 1 CACHE BOOLEAN "Boolean for using PJmalloc library")
ENDIF()

##
## Project declaration
##

project(LMAT VERSION 1.2.4.2018 LANGUAGES C CXX)

#
# Add the right compiler flags
#

# Add -fPIC Wall for convenience.
#CXXFLAGS = -g -O3 -Wall -fopenmp $(INCLUDE) -D$(PJOPT)  \
#            -DIDX_CONFIG=$(SDBIDX) -DTID_SIZE=$(TAXID_SIZE)  -DDBTID_T=$(TID_T) -Wno-deprecated \
#	     -DUSE_SORTED_DB=$(USE_SORTED_DB)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -std=c++0x -DIDX_CONFIG=${SDBIDX} -DTID_SIZE=${TAXID_SIZE} -DDBTID_T=${TID_T} -DUSE_SORTED_DB=${USE_SORTED_DB} -DWITH_PJMALLOC=${PJMALLOC}")
if(USE_BOOST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=${USE_BOOST}")
endif()

# General debug/release compiler flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -g -O0 -Wall")  # -fpermissive")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -Wall")  # -fpermissive")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -Wno-deprecated")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -Wno-deprecated")

# Specific Intel compiler flags
IF(DEFINED IS_INTEL)
  # Enable host architecture optimizations.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xhost")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost")

  # Add further flags depending on build type
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -check=conversions,stack,uninit -ftrapuv -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned -fp-model precise -fp-model source")
  #Requires Intel XE 13.1 
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -ftrapuv -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned -fp-model precise -fp-model source")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -check=conversions,stack,uninit -ftrapuv -fp-model precise -fp-model source")
  
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g -O3 -traceback -debug parallel")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O3 -traceback -debug parallel")
ENDIF()

IF(DEFINED IS_IBMPWR9)
  # Enable host architecture optimizations.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qsmp=omp -qarch=pwr9")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qsmp=omp -qarch=pwr9")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -qsmp=noopt")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -qsmp=noopt")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O5")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O5")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g -O5")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O5")
  set(CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} -qcompact all")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -qcompact all")
ENDIF()

# Specific LLVM/clang compiler flags
if(DEFINED IS_CLANG)
  if(NOT APPLE)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem /software/gcc-4.9.3/include/c++/4.9.3 -I/software/gcc-4.9.3/include/c++/4.9.3/x86_64-unknown-linux-gnu")
  endif()
endif()

#
# Attend external dependencies
#

# Print source and binary dirs
MESSAGE("CMAKE source dir: ${CMAKE_SOURCE_DIR}")
MESSAGE("CMAKE binary dir: ${CMAKE_BINARY_DIR}")

# Set where to look for external packages 
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${CMAKE_MODULE_PATH}) 

MESSAGE("CMAKE_MODULE_PATH source dir: ${CMAKE_MODULE_PATH}")

# Look for OpenMP
FIND_PACKAGE(OpenMP REQUIRED)
IF(OPENMP_FOUND)
  MESSAGE("OpenMP FOUND! We'll try to use it...")
  MESSAGE("  ... with these C flags: ${OpenMP_C_FLAGS}")
  MESSAGE("  ... with these CXX flags: ${OpenMP_CXX_FLAGS}")
  MESSAGE("  ... with these Fortran flags: ${OpenMP_Fortran_FLAGS}")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE(OPENMP_FOUND)
  MESSAGE(FATAL_ERROR "OpenMP not found!")
ENDIF(OPENMP_FOUND)

# Install and check PERM (over JeMalloc)
include(ExternalProject)
ExternalProject_Add(PERM
  PREFIX ${CMAKE_BINARY_DIR}/third-party
  GIT_REPOSITORY https://github.com/khyox/perm-je
  #URL https://computation.llnl.gov/projects/memory-centric-architectures/download/perm-je-0.9.7.tgz
  DOWNLOAD_NO_PROGRESS 1  
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --with-jemalloc-prefix=je CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j 16
  BUILD_IN_SOURCE 1
  LOG_BUILD 0
  INSTALL_COMMAND make install
  LOG_INSTALL 0
)
add_library(jemalloc STATIC IMPORTED)
add_dependencies(jemalloc PERM)
set(PERM_LIB ${CMAKE_BINARY_DIR}/third-party/lib/libjemalloc.a)
set_target_properties(jemalloc PROPERTIES IMPORTED_LOCATION ${PERM_LIB})
message("We'll try to use PERM(JeMalloc)...")
message("  ... with these libraries: ${PERM_LIB}")
# Set include dirs as needed by all_headers.hpp.in inside kmerdb
set(PERM_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/third-party/include)

# Look for thread library
SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads REQUIRED)
IF(THREADS_FOUND)
  MESSAGE("${CMAKE_THREAD_LIBS_INIT} threads lib FOUND! We'll try to use it.")
  IF(CMAKE_USE_PTHREADS_INIT)
    MESSAGE("As preferred, pthreads is to be used.")
  ENDIF(CMAKE_USE_PTHREADS_INIT)
ENDIF(THREADS_FOUND)

# Look for Boost
#FIND_PACKAGE(Boost)
#IF(Boost_FOUND)
#  MESSAGE("Boots ver. ${Boost_VERSION} FOUND! We'll try to use it...")
#  MESSAGE("  ... with these include directories: ${Boost_INCLUDE_DIRS}")
#  MESSAGE("  ... with these libraries: ${Boost_LIBRARIES}")
#ENDIF(Boost_FOUND)

# Look for ZLIB
FIND_PACKAGE(ZLIB REQUIRED)
IF(ZLIB_FOUND)
  MESSAGE("ZLIB ver. ${ZLIB_VERSION_STRING} FOUND! We'll try to use it...")
  MESSAGE("  ... with these include directories: ${ZLIB_INCLUDE_DIRS}")
  MESSAGE("  ... with these libraries: ${ZLIB_LIBRARIES}")
ENDIF(ZLIB_FOUND)

# libgzstream as IMPORTED target library
include(ExternalProject)
ExternalProject_Add(GZ
  PREFIX ${CMAKE_BINARY_DIR}/third-party
  URL http://www.cs.unc.edu/Research/compgeom/gzstream/gzstream.tgz
  DOWNLOAD_NO_PROGRESS 1  
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make CXX=${CMAKE_CXX_COMPILER}
  BUILD_IN_SOURCE 1
  LOG_BUILD 0
  INSTALL_COMMAND mkdir -p ${CMAKE_BINARY_DIR}/third-party/lib
  COMMAND cp <SOURCE_DIR>/libgzstream.a ${CMAKE_BINARY_DIR}/third-party/lib
  COMMAND mkdir -p ${CMAKE_BINARY_DIR}/third-party/include
  COMMAND cp <SOURCE_DIR>/gzstream.h ${CMAKE_BINARY_DIR}/third-party/include
  LOG_INSTALL 0
)
add_library(gzstream STATIC IMPORTED)
add_dependencies(gzstream GZ)
set(GZSTREAM_LIB ${CMAKE_BINARY_DIR}/third-party/lib/libgzstream.a)
set_target_properties(gzstream PROPERTIES IMPORTED_LOCATION ${GZSTREAM_LIB})
message("We'll try to use GZstream...")
message("  ... with these libraries: ${GZSTREAM_LIB}")

# Show summary information before the build
# Print build type and compiler flags
message("<<<>>> SUMMARY <<<>>>")
MESSAGE("The chosen type of build is: ${CMAKE_BUILD_TYPE}")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)
MESSAGE("  >>> with CXX complete flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}}")
MESSAGE("  >>> and linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

# Final message (build successful)
add_custom_target(_LMAT_1_2_4 ALL
                  COMMAND ${CMAKE_COMMAND} -E echo ""
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green "=============================================="
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  " : : :       ··        ··      ·   ·········· "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  " : : :       ···      ···     · ·      ··     "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  " : : :       ·· ··  ·· ··    ·· ··     ··     "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  " : : ······· ··   ··   ··   ·······    ··     "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  " :  ······   ··        ··  ··     ··   ··     "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  "   ·····     ··        ·· ··       ··  ··     "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green "=============================================="
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  "  Livermore  Metagenomics  Analysis  Toolkit  "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green "=============================================="
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan  "  ========== Release 1.2.4 [2018] ==========  "
                  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --green "=============================================="
                  COMMAND ${CMAKE_COMMAND} -E echo ""
                 )

# Dynamically configure a header file with all the headers of the library
set(ALL_HEADERS_HPP ${CMAKE_SOURCE_DIR}/src/kmerdb/all_headers.hpp)
file(GLOB metag_hdr src/kmerdb/*.hpp)
foreach(hpp ${metag_hdr})
  if(NOT ${hpp} STREQUAL ${ALL_HEADERS_HPP})
    set(ALL_HEADERS "${ALL_HEADERS}\n#include \"${hpp}\"")
  endif()
endforeach(hpp)
configure_file(
  "${ALL_HEADERS_HPP}.in"
  "${ALL_HEADERS_HPP}")

#
# Main library and LMAT apps 
#

# Directories that will be searched for headers during compilation
include_directories(${CMAKE_SOURCE_DIR}/include ${PERM_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/kmerdb)

file(GLOB metag_src src/kmerdb/*.cpp)

# From the docs: "Static libraries are always treated as archive targets."
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY src/kmerdb/lib)
# Finally, create it!
add_library(metag STATIC ${metag_src})
set_target_properties(metag PROPERTIES VERSION 1.2.4_2018a SOVERSION 2)
target_link_libraries(metag jemalloc Threads::Threads -ldl gzstream ZLIB::ZLIB ${Boost_LIBRARIES})

# Add the lib as dependency to the final message 
add_dependencies(_LMAT_1_2_4 metag)

# Parse the LMAT apps (add the src programs to the build)
add_subdirectory(src)

# Package-related stuff (use with "cpack")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Livermore Metagenomics Analysis Toolkit")
SET(CPACK_PACKAGE_DESCRIPTION_VENDOR "Lawrence Livermore National Laboratory")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "4")
IF(DEFINED IS_INTEL)
  SET(CPACK_PACKAGE_VERSION_TWEAK "intel")
ELSEIF(DEFINED IS_GNU)
  SET(CPACK_PACKAGE_VERSION_TWEAK "gnu")
ELSEIF(DEFINED IS_CLANG)
  SET(CPACK_PACKAGE_VERSION_TWEAK "clang")
ELSEIF(DEFINED IS_IBMPWR9)
  SET(CPACK_PACKAGE_VERSION_TWEAK "ibmpwr9")
ELSE()
  SET(CPACK_PACKAGE_VERSION_TWEAK "unknown")
ENDIF()
SET(CPACK_GENERATOR "TGZ")

# Get ready for creating a package
INCLUDE(CPack)
