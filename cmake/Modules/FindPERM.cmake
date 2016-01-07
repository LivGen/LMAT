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

find_path(PERM_ROOT_DIR
  NAMES include/jemalloc/pallocator.h
  )

find_library(PERM_LIBRARIES
  NAMES jemalloc
  HINTS ${PERM_ROOT_DIR}/lib
  )

find_path(PERM_INCLUDE_DIRS
  NAMES jemalloc/pallocator.h
  HINTS ${PERM_ROOT_DIR}/include
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PERM DEFAULT_MSG
  PERM_LIBRARIES
  PERM_INCLUDE_DIRS
  )

mark_as_advanced(
  PERM_ROOT_DIR
  PERM_LIBRARIES
  PERM_INCLUDE_DIRS
  )
