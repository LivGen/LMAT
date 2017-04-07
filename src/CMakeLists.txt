# -*- cmake -*-
MESSAGE("CMAKE parsing CMakeLists.txt inside SRC...")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add the test programs to the build
FILE(GLOB lmat_apps ./*.cpp)

# Directories that will be searched for headers during compilation
#INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include ${JEMALLOC_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/kmerdb)

FOREACH(file ${lmat_apps})
  GET_FILENAME_COMPONENT(filename ${file} NAME_WE)
  ADD_EXECUTABLE(${filename} ${file})
  # Add the file as dependency to the final message  
  add_dependencies(_LMAT_1_2_4 ${filename})
  # Link the executable to the lmat library. 
  TARGET_LINK_LIBRARIES(${filename} metag) #${JEMALLOC_LIBRARIES} Threads::Threads -ldl ZLIB::ZLIB ${Boost_LIBRARIES})
ENDFOREACH(file)

