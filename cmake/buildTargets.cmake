file(GLOB_RECURSE vinecopulib_sources src/*.cpp src/*.cc src/*c)
file(GLOB_RECURSE vinecopulib_bicop_headers include/vinecopulib/bicop/*.hpp)
file(GLOB_RECURSE vinecopulib_vinecop_headers include/vinecopulib/vinecop/*.hpp)
file(GLOB_RECURSE vinecopulib_misc_headers include/vinecopulib/misc/*.hpp)
file(GLOB_RECURSE vinecopulib_main_header include/vinecopulib.hpp)
file(GLOB_RECURSE vinecopulib_version_header include/version.hpp)

include_directories(SYSTEM ${external_includes})
include_directories(include)

if (BUILD_SHARED_LIBS)
    add_library(vinecopulib SHARED ${vinecopulib_sources})
else()
    add_library(vinecopulib STATIC ${vinecopulib_sources})
endif()

# TODO: Properly use VINECOPULIB_EXPORT to reduce the number of exports
#include(GenerateExportHeader)
#generate_export_header(vinecopulib)

target_link_libraries(vinecopulib ${external_libs})

set_property(TARGET vinecopulib PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(vinecopulib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS 1)

if(BUILD_TESTING)
    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
    set(unit_tests
            test_all
            test_bicop_parametric
            test_bicop_sanity_checks
            test_bicop_kernel
            test_rvine_matrix
            test_serialization
            test_tools_stats
            test_vinecop_class
            test_vinecop_sanity_checks)

    add_subdirectory(test)
    file(GLOB_RECURSE r_scripts cmake/templates/*R)
    file(COPY ${r_scripts} DESTINATION ${PROJECT_BINARY_DIR}/test)
endif(BUILD_TESTING)

# Related to exports for linux/mac and code coverage
####
# Installation

# Layout. This works for all platforms:
#   * <prefix>/lib/cmake/vinecopulib
#   * <prefix>/lib/
#   * <prefix>/include/
set(config_install_dir "lib/cmake/${PROJECT_NAME}")
set(include_install_dir "include")

set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")

# Configuration
set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")


# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file(
        "${version_config}" COMPATIBILITY SameMajorVersion
)

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * targets_export_name
#   * PROJECT_NAME
configure_package_config_file(
        "cmake/templates/Config.cmake.in"
        "${project_config}"
        INSTALL_DESTINATION "${config_install_dir}"
        PATH_VARS include_install_dir
)

# Targets:
#   * <prefix>/lib/libvinecopulib.dylib
install(
        TARGETS vinecopulib
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        RUNTIME DESTINATION "lib" # on Windows, the dll file is categorised as RUNTIME
)

# Headers:
#   * include/vinecopulib.hpp -> <prefix>/include/vinecopulib.hpp
#   * include/version.hpp -> <prefix>/include/vinecopulib/version.hpp
#   * include/vinecopulib/bicop/*.hpp -> <prefix>/include/vinecopulib/bicop/*.hpp
#   * include/vinecopulib/vinecop/*.hpp -> <prefix>/include/vinecopulib/vinecop/*.hpp
#   * include/vinecopulib/misc/*.hpp -> <prefix>/include/vinecopulib/misc/*.hpp
install(
        FILES ${vinecopulib_main_header}
        DESTINATION "${include_install_dir}"
)
install(
        FILES ${vinecopulib_version_header}
        DESTINATION "${include_install_dir}/vinecopulib"
)
install(
        FILES ${vinecopulib_bicop_headers}
        DESTINATION "${include_install_dir}/vinecopulib/bicop"
)
install(
        FILES ${vinecopulib_vinecop_headers}
        DESTINATION "${include_install_dir}/vinecopulib/vinecop"
)
install(
        FILES ${vinecopulib_misc_headers}
        DESTINATION "${include_install_dir}/vinecopulib/misc"
)

# TODO: Properly use VINECOPULIB_EXPORT to reduce the number of exports
# Export headers:
#   * ${CMAKE_CURRENT_BINARY_DIR}/include_export.h -> <prefix>/include/include_export.h
#install(
#        FILES
#        "${CMAKE_CURRENT_BINARY_DIR}/vinecopulib_export.h"
#        DESTINATION "${include_install_dir}/vinecopulib"
#)

# Config
#   * <prefix>/lib/cmake/vinecopulib/vinecopulibConfig.cmake
#   * <prefix>/lib/cmake/vinecopulib/vinecopulibConfigVersion.cmake
install(
        FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}"
)

# Config
#   * <prefix>/lib/cmake/vinecopulib/vinecopulibTargets.cmake
install(
        EXPORT "${targets_export_name}"
        DESTINATION "${config_install_dir}"
)

# Install the export set for code coverage
if(NOT WIN32 AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING AND CODE_COVERAGE)
    include(cmake/codeCoverage.cmake)
    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/coverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage test_all coverage)
endif()
