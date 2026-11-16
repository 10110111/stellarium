# Locate CFITSIO without downloading dependencies.
find_path(CFITSIO_INCLUDE_DIR fitsio.h)
find_library(CFITSIO_LIBRARY NAMES cfitsio)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CFITSIO REQUIRED_VARS CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)

if(CFITSIO_FOUND AND NOT TARGET CFITSIO::CFITSIO)
    add_library(CFITSIO::CFITSIO UNKNOWN IMPORTED)
    set_target_properties(CFITSIO::CFITSIO PROPERTIES
        IMPORTED_LOCATION "${CFITSIO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INCLUDE_DIR}")
endif()

mark_as_advanced(CFITSIO_INCLUDE_DIR CFITSIO_LIBRARY)
