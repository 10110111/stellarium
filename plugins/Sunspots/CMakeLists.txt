SET(SUNSPOTS_VERSION "0.1.0")
FIND_PACKAGE(CFITSIO QUIET)
IF(NOT CFITSIO_FOUND)
     # CFITSIO calls find_package(ZLIB) and exports its target. Keep the bundled
     # zlib dependency in the build interface so CFITSIO's install export does
     # not refer to Stellarium's source tree or its unexported zlib_stel target.
     IF(NOT USE_SYSTEM_ZLIB)
          SET(ZLIB_INCLUDE_DIR "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/external/zlib>")
          SET(ZLIB_LIBRARY zlib_stel)
          SET(ZLIB_LIBRARIES "$<BUILD_INTERFACE:$<TARGET_FILE:zlib_stel>>")
          SET_TARGET_PROPERTIES(zlib_stel PROPERTIES POSITION_INDEPENDENT_CODE ON)
          IF(NOT TARGET ZLIB::ZLIB)
               ADD_LIBRARY(ZLIB::ZLIB ALIAS zlib_stel)
          ENDIF()
     ENDIF()
     # plugins/CMakeLists.txt defines _REENTRANT for all plugins. CFITSIO uses
     # it to include pthread.h even when its USE_PTHREADS option is off.
     REMOVE_DEFINITIONS(-D_REENTRANT)
     CPMAddPackage(NAME cfitsio
          GIT_REPOSITORY https://github.com/HEASARC/cfitsio.git
          GIT_TAG cfitsio-4.7.0
          OPTIONS
               "TESTS OFF"
               "UTILS OFF"
               "USE_CURL OFF"
               "USE_PTHREADS OFF")
     ADD_DEFINITIONS(-D_REENTRANT)
     IF(NOT USE_SYSTEM_ZLIB)
          # TARGET_FILE in a link expression does not add a build dependency.
          ADD_DEPENDENCIES(cfitsio zlib_stel)
     ENDIF()
     ADD_LIBRARY(CFITSIO::CFITSIO ALIAS cfitsio)
ENDIF()
ADD_DEFINITIONS(-DSUNSPOTS_PLUGIN_VERSION="${SUNSPOTS_VERSION}")
ADD_DEFINITIONS(-DSUNSPOTS_PLUGIN_LICENSE="GNU GPLv2 or later")

ADD_SUBDIRECTORY( src )
