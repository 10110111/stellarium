########### install files ###############

INSTALL(DIRECTORY ./ DESTINATION share/${PACKAGE}/textures
        FILES_MATCHING PATTERN "*.png" PATTERN "atmosphere/*.dat"
        PATTERN "CMakeFiles" EXCLUDE )


