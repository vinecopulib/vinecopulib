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
  find_package(Eigen3 REQUIRED)
  if (EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR} (found suitable version \"${Eigen3_VERSION}\")")
  else()
    message(FATAL_ERROR "Could not find Eigen3")
  endif()
endif()

# Check if Boost_INCLUDE_DIRS is defined and if not, try to find it
if(NOT DEFINED Boost_INCLUDE_DIRS)
  # try to find Boost in CONFIG mode first
  find_package(Boost 1.56 CONFIG)
  if (Boost_FOUND)
    message(STATUS "Found Boost: ${Boost_DIR} (found suitable version \"${Boost_VERSION}\")")  
  else ()
    # fallback to MODULE mode
    find_package(Boost 1.56 MODULE REQUIRED)
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
      GIT_TAG        c837460 
    )
    FetchContent_MakeAvailable(wdm)
    set(wdm_INCLUDE_DIRS "${wdm_SOURCE_DIR}/include")
  else()
    message(STATUS "Found wdm: ${wdm_INCLUDE_DIRS} (found suitable version \"${wdm_VERSION}\")")
  endif()
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

# Set all the external dependencies
set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})

if(BUILD_DOC)
  # Find doxygen and configure if found
  find_package(Doxygen REQUIRED)
  configure_file(
          ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
      )
  add_custom_target(doc
          ${DOXYGEN_EXECUTABLE}
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif()
