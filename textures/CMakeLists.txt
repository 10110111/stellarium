########### install files ###############

INSTALL(DIRECTORY ./ DESTINATION ${SDATALOC}/textures
        FILES_MATCHING PATTERN "*.png" PATTERN "atmosphere/*.dat"
        PATTERN "CMakeFiles" EXCLUDE )


