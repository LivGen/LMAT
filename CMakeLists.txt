# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

#
# Compiler family selection
#

# Choose compiler family by env var COMPILER_FAMILY
IF(COMPILER_FAMILY STREQUAL "intel")
  SET(IS_INTEL TRUE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by INTEL Compilers...")
ELSE()
  SET(IS_INTEL FALSE)
  MESSAGE("CMAKE parsing CMakeLists.txt for LMAT by GNU Compilers...")
ENDIF()

IF(IS_INTEL)
  # Force Intel Compilers (do it before PROJECT!)
  INCLUDE(CMakeForceCompiler)
  MESSAGE("Enforcing Intel compilers!")
  set(CMAKE_C_COMPILER icc)
  set(CMAKE_CXX_COMPILER icpc)
  set(CMAKE_Fortran_COMPILER ifort)
ELSE()
  set(CMAKE_C_COMPILER gcc)
  set(CMAKE_CXX_COMPILER g++)
  set(CMAKE_Fortran_COMPILER gfortran)  
ENDIF(IS_INTEL)

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
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
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


##
## Project declaration
##

project(LMAT VERSION 1.2.4.0 LANGUAGES C CXX)


#
# Add the right compiler flags
#

# Add -fPIC Wall for convenience.
#CXXFLAGS = -g -O3 -Wall -fopenmp $(INCLUDE) -D$(PJOPT)  \
#            -DIDX_CONFIG=$(SDBIDX) -DTID_SIZE=$(TAXID_SIZE)  -DDBTID_T=$(TID_T) -Wno-deprecated \
#	     -DUSE_SORTED_DB=$(USE_SORTED_DB)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -std=c++0x -DIDX_CONFIG=${SDBIDX} -DTID_SIZE=${TAXID_SIZE} -DDBTID_T=${TID_T} -DUSE_SORTED_DB=${USE_SORTED_DB}")
if(USE_BOOST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_BOOST=${USE_BOOST}")
endif()

# General debug/release compiler flags
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -g -O0 -Wall")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g -O0 -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -Wno-deprecated")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -Wno-deprecated")

# Specific Intel compiler flags
IF(IS_INTEL)
  # Enable host architecture optimizations.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xhost")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost")

  # Add further flags depending on build type
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -check=conversions,stack,uninit -ftrapuv -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned -fp-model precise -fp-model source")
  #Requires Intel XE 13.1 
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -check-uninit -ftrapuv -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned -fp-model precise -fp-model source")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -traceback -rdynamic -debug parallel -check=conversions,stack,uninit -ftrapuv -fp-model precise -fp-model source")
  
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -g -O3 -traceback -debug parallel")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O3 -traceback -debug parallel")
ENDIF()


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
ENDIF(OPENMP_FOUND)

# Install and check PERM (over JeMalloc)
include(ExternalProject)
ExternalProject_Add(PERM
  PREFIX ${CMAKE_BINARY_DIR}/third-party
  URL http://computation.llnl.gov/project/perma/downloads/perm-je-0.9.7.tgz
  DOWNLOAD_NO_PROGRESS 1  
  CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> CXX=${CMAKE_CXX_COMPILER}
  BUILD_COMMAND make -j 16
  BUILD_IN_SOURCE 1
  LOG_BUILD 0
  INSTALL_COMMAND make install
  LOG_INSTALL 0
)
add_library(jemalloc STATIC IMPORTED)
add_dependencies(jemalloc PERM)
#target_link_libraries(LMAT jemalloc)
set(PERM_ROOT_DIR ${CMAKE_BINARY_DIR}/third-party)
set(PERM_LIBRARIES ${PERM_ROOT_DIR}/lib/libjemalloc.a)
set(PERM_INCLUDE_DIRS ${PERM_ROOT_DIR}/include)
set(PERM_FOUND 1)
mark_as_advanced(FORCE
  PERM_ROOT_DIR
  PERM_LIBRARIES
  PERM_INCLUDE_DIRS
  PERM_FOUND
  )
#FIND_PACKAGE(PERM REQUIRED)
IF(PERM_FOUND)
  MESSAGE("We'll try to use PERM(JeMalloc)...")
  MESSAGE("  ... from this root dir: ${PERM_ROOT_DIR}")
  MESSAGE("  ... with these include directories: ${PERM_INCLUDE_DIRS}")
  MESSAGE("  ... with these libraries: ${PERM_LIBRARIES}")
ENDIF(PERM_FOUND)

# Look for thread library
SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads REQUIRED)
IF(THREADS_FOUND)
  MESSAGE("${CMAKE_THREAD_LIBS_INIT} threads lib FOUND! We'll try to use it.")
  IF(CMAKE_USE_PTHREADS_INIT)
    MESSAGE("As preferred, pthreads is to be used.")
  ENDIF(CMAKE_USE_PTHREADS_INIT)
ENDIF(THREADS_FOUND)
##### ADD THIS
##### target_link_libraries(my_app Threads::Threads)

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
##### ADD THIS
##### target_link_libraries(my_app ZLIB::ZLIB)

# libgzstream as IMPORTED target library
#add_library(gzstream STATIC IMPORTED)
#set_property(TARGET gzstream PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/libgzstream.a)
# Install and check PERM (over JeMalloc)
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
#target_link_libraries(LMAT jemalloc)
set(GZ_ROOT_DIR ${CMAKE_BINARY_DIR}/third-party)
set(GZ_LIBRARIES ${GZ_ROOT_DIR}/lib/libgzstream.a)
set(GZ_INCLUDE_DIRS ${GZ_ROOT_DIR}/include)
set(GZ_FOUND 1)
mark_as_advanced(FORCE
  GZ_ROOT_DIR
  GZ_LIBRARIES
  GZ_INCLUDE_DIRS
  GZ_FOUND
  )
#FIND_PACKAGE(PERM REQUIRED)
IF(GZ_FOUND)
  MESSAGE("We'll try to use gzstream...")
  MESSAGE("  ... from this root dir: ${GZ_ROOT_DIR}")
  MESSAGE("  ... with these include directories: ${GZ_INCLUDE_DIRS}")
  MESSAGE("  ... with these libraries: ${GZ_LIBRARIES}")
ENDIF(GZ_FOUND)


#
# Show summary information before the build
#

# Print build type and compiler flags
MESSAGE("The chosen type of build is: ${CMAKE_BUILD_TYPE}")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPERCASE)
MESSAGE("  with CXX complete flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPERCASE}}")
MESSAGE("  and linker flags: ${CMAKE_EXE_LINKER_FLAGS}")


#
# Final message (build successful)
#
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
                  COMMAND ${CMAKE_COMMAND} -E echo ""
                 )


#
# Create metagen lib in the right place
#

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

# Directories that will be searched for headers during compilation
include_directories(${CMAKE_SOURCE_DIR}/include ${PERM_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GZ_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/kmerdb)

file(GLOB metag_src src/kmerdb/*.cpp)

# From the docs: "Static libraries are always treated as archive targets."
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY src/kmerdb/lib)
# Finally, create it!
add_library(metag STATIC ${metag_src})
set_target_properties(metag PROPERTIES VERSION 1.2.4 SOVERSION 2)
target_link_libraries(metag ${PERM_LIBRARIES} Threads::Threads -ldl ${GZ_LIBRARIES} ZLIB::ZLIB ${Boost_LIBRARIES})
# Add the lib as dependency to the final message 
add_dependencies(_LMAT_1_2_4 metag)


#
# Parse the LMAT apps
#

# Add the src programs to the build
add_subdirectory(src)


#
# Package-related stuff (use with "cpack")
#

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Livermore Metagenomics Analysis Toolkit")
SET(CPACK_PACKAGE_DESCRIPTION_VENDOR "Lawrence Livermore National Laboratory")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "4")
IF(IS_INTEL)
  SET(CPACK_PACKAGE_VERSION_TWEAK "intel")
ELSE()
  SET(CPACK_PACKAGE_VERSION_TWEAK "gnu")
ENDIF()
SET(CPACK_GENERATOR "TGZ")

# Get ready for creating a package
INCLUDE(CPack)
