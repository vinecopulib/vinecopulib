include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

add_subdirectory(src)
if(BUILD_TESTING)
    set(unit_tests
            test_all
            test_bicop_parametric
            test_bicop_class
            test_bicop_select
            test_bicop_trafokernel
            test_rvine_matrix
            test_vinecop_class)
    add_subdirectory(test)
endif(BUILD_TESTING)

# Add all targets to the build-tree export set
export(TARGETS  vinecopulib FILE "${PROJECT_BINARY_DIR}/vinecopulibTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE vinecopulib)

# Create the vinecopulib-config.cmake and vinecopulib-config-version.cmake files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(cmake/templates/vinecopulib-config.cmake.in
        "${PROJECT_BINARY_DIR}/vinecopulib-config.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${VINECOPULIB_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/templates/vinecopulib-config.cmake.in
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vinecopulib-config.cmake" @ONLY)
# ... for both
configure_file(cmake/templates/vinecopulib-config-version.cmake.in
        "${PROJECT_BINARY_DIR}/vinecopulib-config-version.cmake" @ONLY)

# Install the vinecopulib-config.cmake and vinecopulib-config-version.cmake
install(FILES
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vinecopulib-config.cmake"
        "${PROJECT_BINARY_DIR}/vinecopulib-config-version.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}")

# Install the export set for use with the install-tree
install(EXPORT vinecopulibTargets DESTINATION "${INSTALL_CMAKE_DIR}")

# Install the export set for code coverage
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING)
    include(cmake/codeCoverage.cmake)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/coverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage test_all coverage)
endif()