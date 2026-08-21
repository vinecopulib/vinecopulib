if (POLICY CMP0074)
  # find_package() uses <PackageName>_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif ()

if (POLICY CMP0144)
  # find_package() uses upper-case <PACKAGENAME>_ROOT variables.
  cmake_policy(SET CMP0144 NEW)
endif ()

# Find the main dependencies

# Check if EIGEN3_INCLUDE_DIR is defined and if not, try to find it
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 REQUIRED CONFIG)
  if (Eigen3_FOUND)
    message(STATUS "Found Eigen3: ${Eigen3_DIR} (found suitable version \"${Eigen3_VERSION}\")")
  else()
    message(FATAL_ERROR "Could not find Eigen3")
  endif()
  # Eigen 5.x exposes its include path only through the target and no longer
  # sets EIGEN3_INCLUDE_DIR, which external_includes below needs.
  if(TARGET Eigen3::Eigen)
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen
                        INTERFACE_INCLUDE_DIRECTORIES)
  endif()
endif()

# Check if Boost_INCLUDE_DIRS is defined and if not, try to find it
if(NOT DEFINED Boost_INCLUDE_DIRS)
  # try to find Boost in CONFIG mode first
  find_package(Boost 1.75 CONFIG)
  if (Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_DIR} (found suitable version \"${Boost_VERSION}\")")  
  else ()
    # fallback to MODULE mode
    find_package(Boost 1.75 MODULE REQUIRED)
  endif ()
endif()

find_package(Threads                      REQUIRED)

# Check if wdm_INCLUDE_DIRS is defined and if not, try to find it
if(NOT DEFINED wdm_INCLUDE_DIRS)
  # Download if not found
  find_package(wdm 0.2.6 QUIET)
  if(NOT wdm_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      wdm
      GIT_REPOSITORY https://github.com/tnagler/wdm.git
      GIT_TAG        6b85eebc3e0097672be1047f6f43848bc77acb52 # tnagler/wdm#17
    )
    FetchContent_MakeAvailable(wdm)
    set(wdm_INCLUDE_DIRS "${wdm_SOURCE_DIR}/include")
  else()
    message(STATUS "Found wdm: ${wdm_INCLUDE_DIRS} (found suitable version \"${wdm_VERSION}\")")
  endif()
endif()

# Chatterjee's xi carries wdm's 0.2.6 version number, so neither the version
# check above nor a caller-supplied wdm_INCLUDE_DIRS can rule out a wdm that
# predates it.
if(NOT EXISTS "${wdm_INCLUDE_DIRS}/wdm/cxi.hpp")
  message(FATAL_ERROR "The wdm at ${wdm_INCLUDE_DIRS} is missing wdm/cxi.hpp, "
          "which the \"cxi\" tree criterion needs. Point wdm_DIR or "
          "wdm_INCLUDE_DIRS at wdm 6b85eeb or later, or unset both to have it "
          "fetched at configure time.")
endif()

# The build links wdm by target name, which find_package and FetchContent both
# define but a caller-supplied wdm_INCLUDE_DIRS does not. Without this the name
# is taken for a library to search for, and every link fails on -lwdm.
if(NOT TARGET wdm)
  add_library(wdm INTERFACE IMPORTED GLOBAL)
  set_target_properties(wdm PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${wdm_INCLUDE_DIRS}")
endif()

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
