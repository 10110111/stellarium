INCLUDE_DIRECTORIES(
     .
     gui
     ${CMAKE_BINARY_DIR}/plugins/Sunspots/src
     ${CMAKE_BINARY_DIR}/plugins/Sunspots/src/gui
)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/src)

SET(Sunspots_SRCS
     Sunspots.hpp
     Sunspots.cpp
     SunspotsImage.hpp
     SunspotsImage.cpp
     gui/SunspotsDialog.hpp
     gui/SunspotsDialog.cpp
)

SET(Sunspots_UIS
     gui/sunspotsDialog.ui
)

SET(Sunspots_RES
     ../resources/Sunspots.qrc
)

IF (${QT_VERSION_MAJOR} EQUAL "5")
     QT5_WRAP_UI(Sunspots_UIS_H ${Sunspots_UIS})
     QT5_ADD_RESOURCES(Sunspots_RES_CXX ${Sunspots_RES})
ELSE()
     QT_WRAP_UI(Sunspots_UIS_H ${Sunspots_UIS})
     QT_ADD_RESOURCES(Sunspots_RES_CXX ${Sunspots_RES})
ENDIF()

ADD_LIBRARY(Sunspots-static STATIC ${Sunspots_SRCS} ${Sunspots_RES_CXX} ${Sunspots_UIS_H})
TARGET_INCLUDE_DIRECTORIES(Sunspots-static PRIVATE ${OPENJPEG_INCLUDE_DIRS})
# Windows static OpenJPEG consumers must not import its symbols from a DLL.
GET_TARGET_PROPERTY(SUNSPOTS_OPENJPEG_TYPE OpenJPEG::OpenJPEG TYPE)
IF(SUNSPOTS_OPENJPEG_TYPE STREQUAL "STATIC_LIBRARY")
     TARGET_COMPILE_DEFINITIONS(Sunspots-static PRIVATE OPJ_STATIC)
ENDIF()
TARGET_LINK_LIBRARIES(Sunspots-static Qt${QT_VERSION_MAJOR}::Core Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::OpenGL Qt${QT_VERSION_MAJOR}::Network CFITSIO::CFITSIO OpenJPEG::OpenJPEG)
SET_TARGET_PROPERTIES(Sunspots-static PROPERTIES OUTPUT_NAME "Sunspots")
SET_TARGET_PROPERTIES(Sunspots-static PROPERTIES COMPILE_FLAGS "-DQT_STATICPLUGIN")

ADD_DEPENDENCIES(AllStaticPlugins Sunspots-static)

SET_TARGET_PROPERTIES(Sunspots-static PROPERTIES FOLDER "plugins/Sunspots")
