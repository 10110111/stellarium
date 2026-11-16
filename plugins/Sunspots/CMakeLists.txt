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
          VERSION 4.7.0
          URL https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.7.0.tar.gz
          URL_HASH SHA256=ce573bbea8e75b429f8c3d3e86498741ba3dc9628a1530d2f65268397ad059e8
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
# Qt's image format plugins do not consistently provide JPEG2000 support.
FIND_PACKAGE(OpenJPEG CONFIG QUIET)
IF(NOT OpenJPEG_FOUND)
     # OpenJPEG inserts these legacy output paths into the global CMake cache.
     # Do not let its defaults change where Stellarium and subsequent plugins are built.
     SET(SUNSPOTS_OPENJPEG_NEW_OUTPUT_PATHS)
     FOREACH(OUTPUT_PATH EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)
          IF(NOT DEFINED ${OUTPUT_PATH})
               LIST(APPEND SUNSPOTS_OPENJPEG_NEW_OUTPUT_PATHS ${OUTPUT_PATH})
          ENDIF()
     ENDFOREACH()
     CPMAddPackage(NAME openjpeg
          VERSION 2.5.3
          URL https://github.com/uclouvain/openjpeg/archive/refs/tags/v2.5.3.tar.gz
          URL_HASH SHA256=368fe0468228e767433c9ebdea82ad9d801a3ad1e4234421f352c8b06e7aa707
          OPTIONS
               "BUILD_CODEC OFF"
               "BUILD_SHARED_LIBS OFF"
               "BUILD_STATIC_LIBS ON"
               "BUILD_DOC OFF"
               "BUILD_TESTING OFF"
               "BUILD_UNIT_TESTS OFF"
               "BUILD_JPIP OFF"
               "BUILD_VIEWER OFF"
               "BUILD_JAVA OFF")
     FOREACH(OUTPUT_PATH ${SUNSPOTS_OPENJPEG_NEW_OUTPUT_PATHS})
          UNSET(${OUTPUT_PATH} CACHE)
     ENDFOREACH()
     SET(OPENJPEG_INCLUDE_DIRS
          ${openjpeg_SOURCE_DIR}/src/lib/openjp2
          ${openjpeg_BINARY_DIR}/src/lib/openjp2)
     SET_TARGET_PROPERTIES(openjp2 PROPERTIES POSITION_INDEPENDENT_CODE ON AUTOMOC OFF)
ENDIF()
IF(NOT TARGET OpenJPEG::OpenJPEG)
     ADD_LIBRARY(OpenJPEG::OpenJPEG ALIAS openjp2)
ENDIF()

ADD_DEFINITIONS(-DSUNSPOTS_PLUGIN_VERSION="${SUNSPOTS_VERSION}")
ADD_DEFINITIONS(-DSUNSPOTS_PLUGIN_LICENSE="GNU GPLv2 or later")

ADD_SUBDIRECTORY( src )
