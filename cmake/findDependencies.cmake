# The build links these by target name, so a caller-supplied path has to define
# the target too. Call after the find_package that may already define it.
function(vinecopulib_add_header_only_target target include_dirs_var)
  if(TARGET ${target})
    return()
  endif()
  if(NOT ${include_dirs_var})
    message(FATAL_ERROR
            "${target} is not defined and ${include_dirs_var} is empty. Set "
            "${include_dirs_var} to the directory holding the headers, or "
            "leave it unset and let find_package look for the package.")
  endif()
  add_library(${target} INTERFACE IMPORTED GLOBAL)
  set_target_properties(${target} PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${${include_dirs_var}}")
endfunction()

# Find the main dependencies

# Check if EIGEN3_INCLUDE_DIR is defined and if not, try to find it
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 REQUIRED CONFIG)
  message(STATUS "Found Eigen3: ${Eigen3_DIR} (found suitable version \"${Eigen3_VERSION}\")")
  # Eigen 5.x sets no EIGEN3_INCLUDE_DIR; external_includes below needs one.
  get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen
                      INTERFACE_INCLUDE_DIRECTORIES)
endif()
vinecopulib_add_header_only_target(Eigen3::Eigen EIGEN3_INCLUDE_DIR)

# Check if Boost_INCLUDE_DIRS is defined and if not, try to find it
if(NOT DEFINED Boost_INCLUDE_DIRS)
  # CONFIG only: CMP0167 removed the FindBoost module, and Boost has shipped a
  # config package since 1.70, below the 1.75 required here.
  find_package(Boost 1.75 REQUIRED CONFIG)
  message(STATUS "Found Boost: ${Boost_DIR} (found suitable version \"${Boost_VERSION}\")")
endif()
vinecopulib_add_header_only_target(Boost::headers Boost_INCLUDE_DIRS)

find_package(Threads                      REQUIRED)

# Check if wdm_INCLUDE_DIRS is defined and if not, try to find it
if(NOT DEFINED wdm_INCLUDE_DIRS)
  # Download if not found
  # 0.3.0 for Chatterjee's xi; do not lower.
  find_package(wdm 0.3.0 QUIET)
  if(NOT wdm_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      wdm
      GIT_REPOSITORY https://github.com/tnagler/wdm.git
      GIT_TAG        v0.3.0
    )
    FetchContent_MakeAvailable(wdm)
    set(wdm_INCLUDE_DIRS "${wdm_SOURCE_DIR}/include")
  else()
    message(STATUS "Found wdm: ${wdm_INCLUDE_DIRS} (found suitable version \"${wdm_VERSION}\")")
  endif()
endif()

vinecopulib_add_header_only_target(wdm wdm_INCLUDE_DIRS)

# Ensure R is available and download googlestest
if(BUILD_TESTING)
  include(cmake/findR.cmake                 REQUIRED)

  # Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  include(FetchContent)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        f8d7d77c06936315286eb55f8de22cd23c188571 # release-1.14.0
  )

  # Download and configure googletest
  FetchContent_MakeAvailable(googletest)
endif()

# Download google/benchmark for the benchmark suite
if(VINECOPULIB_BUILD_BENCHMARKS)
  include(FetchContent)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
  set(BENCHMARK_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        96afad55c79e02f5dfca1374e772c2be72ba631b # v1.9.1
  )
  FetchContent_MakeAvailable(googlebenchmark)
endif()

# Set all the external dependencies
set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})

if(VINECOPULIB_BUILD_DOC)
  # Find doxygen and configure if found
  find_package(Doxygen REQUIRED)
  configure_file(
          ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
      )
  # The m.css variant only overrides a few tags and @INCLUDEs the above.
  configure_file(
          ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile-mcss.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile-mcss @ONLY
      )
  add_custom_target(doc
          ${DOXYGEN_EXECUTABLE}
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
  # The snippets the pages reference must compile before the pages are built.
  add_dependencies(doc doc_snippets)
endif()
