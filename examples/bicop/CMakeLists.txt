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
find_package(Threads                      REQUIRED)

# Set required variables for includes and libraries
# In the second line
#   * VINECOPULIB_LIBRARIES is needed if vinecopulib has been built as a
#     shared lib (does nothing otherwise).
#   * CMAKE_THREAD_LIBS_INIT is needed for some linux systems
#     (but does nothing on OSX/Windows).
set(external_includes ${VINECOPULIB_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(external_libs ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Include subdirectory with project sources
add_subdirectory(src)
