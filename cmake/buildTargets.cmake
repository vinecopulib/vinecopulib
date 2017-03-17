include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

file(GLOB_RECURSE vinecopulib_sources src/*.cpp src/*.cc src/*c)
file(GLOB_RECURSE vinecopulib_headers include/*.h include/*.hpp)

include_directories(${external_includes} include)

if (BUILD_SHARED_LIBS)
    add_library(vinecopulib SHARED ${vinecopulib_sources})
else()
    add_library(vinecopulib STATIC ${vinecopulib_sources})
endif()

include(GenerateExportHeader)
generate_export_header(vinecopulib)

target_link_libraries(vinecopulib ${external_libs})

set_property(TARGET vinecopulib PROPERTY POSITION_INDEPENDENT_CODE ON)

set_target_properties(vinecopulib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS 1)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

if(BUILD_TESTING)
    set(unit_tests
            test_all
            test_bicop_parametric
            test_bicop_sanity_checks
            test_bicop_select
            test_bicop_trafokernel
            test_rvine_matrix
            test_vinecop_class)
    add_subdirectory(test)
endif(BUILD_TESTING)

# Related to code coverage and exports in linux
if (NOT WIN32)
    ####
    # Installation

    # Layout. This works for all platforms:
    #   * <prefix>/lib/cmake/<PROJECT-NAME>
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
    #   * header location after install: <prefix>/include/vinecopulib/*.hpp
    #   * headers can be included by C++ code `#include <vinecopulib/*.hpp>`
    install(
            TARGETS vinecopulib
            EXPORT "${targets_export_name}"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib"
            RUNTIME DESTINATION "bin"
    )

    # Headers:
    #   * include/vinecopulib/*.hpp -> <prefix>/include/vinecopulib/*.hpp
    install(
            FILES ${vinecopulib_headers}
            DESTINATION "${include_install_dir}/vinecopulib"
    )

    # Export headers:
    #   * ${CMAKE_CURRENT_BINARY_DIR}/include_export.h -> <prefix>/include/include_export.h
    install(
            FILES
            "${CMAKE_CURRENT_BINARY_DIR}/vinecopulib_export.h"
            DESTINATION "${include_install_dir}/vinecopulib"
    )

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
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND BUILD_TESTING)
        include(cmake/codeCoverage.cmake)
        file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/coverage)
        setup_target_for_coverage(${PROJECT_NAME}_coverage test_all coverage)
    endif()
endif()
