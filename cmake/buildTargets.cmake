if (VINECOPULIB_SHARED_LIB)
    add_library(vinecopulib ${vinecopulib_sources})
    target_include_directories(vinecopulib PRIVATE ${vinecopulib_includes})
    set_property(TARGET vinecopulib PROPERTY POSITION_INDEPENDENT_CODE ON)
    set_target_properties(vinecopulib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS 1)
    target_link_libraries(vinecopulib PUBLIC Eigen3::Eigen wdm Boost::boost ${CMAKE_THREAD_LIBS_INIT})
    # non windows
    if (NOT WIN32)
        target_compile_options(vinecopulib PRIVATE -Wno-maybe-uninitialized) # Boost triggers this warning in strict mode
    endif()
else()
    add_library(vinecopulib INTERFACE)
    target_link_libraries(vinecopulib INTERFACE Eigen3::Eigen wdm Boost::boost ${CMAKE_THREAD_LIBS_INIT})
endif()

target_include_directories(vinecopulib INTERFACE $<BUILD_INTERFACE:${vinecopulib_includes}>)
target_include_directories (vinecopulib INTERFACE $<INSTALL_INTERFACE:include>)

if(BUILD_TESTING)

    set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
    set(unit_tests
            test_all
            test_bicop_parametric
            test_bicop_sanity_checks
            test_bicop_kernel
            test_bicop_select
            test_rvine_structure
            test_tools_bobyqa
            test_tools_stats
            test_vinecop_class
            test_vinecop_sanity_checks
            test_weights
            test_discrete)

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
if (VINECOPULIB_SHARED_LIB)
    install(TARGETS vinecopulib
            EXPORT "${targets_export_name}"
            LIBRARY DESTINATION "lib"
            ARCHIVE DESTINATION "lib"
            RUNTIME DESTINATION "bin" # on Windows, the dll file is categorised as RUNTIME
    )
else()
    install(TARGETS vinecopulib EXPORT "${targets_export_name}")
endif()


# Headers:
#   * include/vinecopulib.hpp -> <prefix>/include/vinecopulib.hpp
#   * include/version.hpp -> <prefix>/include/vinecopulib/version.hpp
#   * include/vinecopulib/bicop/*.hpp -> <prefix>/include/vinecopulib/bicop/*.hpp
#   * include/vinecopulib/vinecop/*.hpp -> <prefix>/include/vinecopulib/vinecop/*.hpp
#   * include/vinecopulib/misc/*.hpp -> <prefix>/include/vinecopulib/misc/*.hpp
#   * wdm/include/*.hpp -> <prefix>/include/vinecopulib/wdm/*.hpp
install(
        FILES ${main_hpp}
        DESTINATION "${include_install_dir}"
)
install(
        FILES ${version_hpp}
        DESTINATION "${include_install_dir}/vinecopulib"
)
install(
        FILES ${bicop_hpp}
        DESTINATION "${include_install_dir}/vinecopulib/bicop"
)
install(
        FILES ${vinecop_hpp}
        DESTINATION "${include_install_dir}/vinecopulib/vinecop"
)
install(
        FILES ${misc_hpp}
        DESTINATION "${include_install_dir}/vinecopulib/misc"
)
if (NOT VINECOPULIB_SHARED_LIB)
    install(
            FILES ${vinecopulib_bicop_ipp}
            DESTINATION "${include_install_dir}/vinecopulib/bicop/implementation"
    )
    install(
            FILES ${vinecopulib_vinecop_ipp}
            DESTINATION "${include_install_dir}/vinecopulib/vinecop/implementation"
    )
    install(
            FILES ${vinecopulib_misc_ipp}
            DESTINATION "${include_install_dir}/vinecopulib/misc/implementation"
    )
endif()

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
