cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project (Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
if (NOT WIN32)
  set(CMAKE_CXX_FLAGS "-std=gnu++11 -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O2 -DNDEBUG")
endif()

# Find vinecopulib package and dependencies
find_package(vinecopulib                  REQUIRED)
find_package(Boost 1.56                   REQUIRED)
include(cmake/findEigen3.cmake            REQUIRED)
find_package(NLopt)
if ( NLOPT_INCLUDE_DIRS )
  # this is needed for GCC 6 because we later use these as SYSTEM headers
  # see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=70129
  #     https://issues.apache.org/jira/browse/THRIFT-3828
  #     https://gitlab.kitware.com/cmake/cmake/issues/16291
  get_filename_component(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS} REALPATH)
else ( NLOPT_INCLUDE_DIRS )
  # for NLopt versions installed with apt-get/brew
  message(STATUS "NLopt not found via find_package, try locating it in standard repos.")
  include(cmake/findNLopt.cmake            REQUIRED)
endif()

# Set required variables for includes and libraries
set(external_includes ${VINECOPULIB_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(external_libs ${VINECOPULIB_LIBRARIES} ${NLOPT_LIBRARIES} ${Boost_LIBRARIES})

# Include subdirectory with project sources
add_subdirectory(src)
