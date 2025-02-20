set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(GettextPo REQUIRED)

add_executable(skyculture-converter
                Utils.cpp
                NamesOldLoader.cpp
                AsterismOldLoader.cpp
                DescriptionOldLoader.cpp
                ConstellationOldLoader.cpp
                main.cpp
              )
target_link_libraries(skyculture-converter Qt::Core Qt::Gui GettextPo::GettextPo)

if(NOT APPLE) # FIXME: what to do on macOS?
     install(TARGETS skyculture-converter DESTINATION bin)
endif()
