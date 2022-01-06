# Find the main dependencies
find_package(Eigen3                       REQUIRED)
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        703bd9caab50b139428cea1aaff9974ebee5742e # release-1.10.0
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_Declare(
  wdm
  URL https://github.com/tnagler/wdm/archive/refs/heads/dev.zip
)

# Find wdm and download if not found
# Download googlestest
find_package(wdm QUIET)
if(NOT wdm_FOUND)
  if(BUILD_TESTING)
    FetchContent_MakeAvailable(wdm googletest)
  else()
    FetchContent_MakeAvailable(wdm)
  endif()
  set(wdm_INCLUDE_DIRS "${wdm_SOURCE_DIR}/include")
elseif(BUILD_TESTING)
  FetchContent(googletest)
endif()

# Set all the external dependencies
set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})

# Find doxygen and configure if found
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
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
endif (DOXYGEN_FOUND)
