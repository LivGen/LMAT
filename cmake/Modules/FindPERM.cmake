# - Try to find PERM (over jemalloc) headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(PERM)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  PERM_ROOT_DIR Set this variable to the root installation of
#                    PERM if the module has problems finding
#                    the proper installation path.
#
# Variables defined by this module:
#
#  PERM_FOUND             System has PERM libs/headers
#  PERM_LIBRARIES         The PERM library/libraries
#  PERM_INCLUDE_DIRS      The location of PERM headers

set(PERM_ROOT_DIR ${CMAKE_BINARY_DIR}/third-party)

#find_path(PERM_ROOT_DIR
#  NAMES include/jemalloc/pallocator.h
#  )

find_library(PERM_LIBRARY
  NAMES jemalloc perm-je
  PATHS ${PERM_ROOT_DIR}/lib
  )

find_path(PERM_INCLUDE_DIR
  NAMES jemalloc/pallocator.h
  PATH ${PERM_ROOT_DIR}/include
  )

MESSAGE("  ... from this root dir: ${PERM_ROOT_DIR}")
MESSAGE("  ... with these include directories: ${PERM_INCLUDE_DIRS}")
MESSAGE("  ... with these libraries: ${PERM_LIBRARIES}")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set <PackageName>_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(PERM DEFAULT_MSG PERM_INCLUDE_DIR PERM_LIBRARY)

if (PERM_FOUND)
    set(PERM_LIBRARIES ${PERM_LIBRARY} )
    set(PERM_INCLUDE_DIRS ${PERM_INCLUDE_DIR} )
    set(PERM_DEFINITIONS )
endif()

mark_as_advanced(
  PERM_ROOT_DIR
  PERM_LIBRARIES
  PERM_INCLUDE_DIRS
  )
