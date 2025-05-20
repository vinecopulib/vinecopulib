cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)

project (Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
if (NOT WIN32)
  set(CMAKE_CXX_FLAGS "-Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O3 -march=native -DNDEBUG")
endif()

if (POLICY CMP0074)
  # find_package() uses <PackageName>_ROOT variables
  cmake_policy(SET CMP0074 NEW)
endif ()

if (POLICY CMP0144)
  # find_package() uses upper-case <PACKAGENAME>_ROOT variables.
  cmake_policy(SET CMP0144 NEW)
endif ()

# Find vinecopulib package and dependencies
find_package(vinecopulib REQUIRED)

# Include subdirectory with project sources
add_subdirectory(src)
