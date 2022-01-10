# Find the main dependencies
find_package(Eigen3                       REQUIRED)
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        6b74da4757a549563d7c37c8fae3e704662a043b # release-1.10.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_Declare(
  wdm
  GIT_REPOSITORY https://github.com/tnagler/wdm.git
  GIT_TAG        bf49a5709b795c4501c25a99ecb18e838cbab631
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
