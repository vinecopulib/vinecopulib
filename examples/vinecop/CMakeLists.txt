cmake_minimum_required(VERSION 3.20...4.0)

set(CMAKE_CXX_STANDARD 17)

project (Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
if (NOT WIN32)
    set(CMAKE_CXX_FLAGS "-Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O3 -DNDEBUG")
endif()

# Find vinecopulib package and dependencies
find_package(vinecopulib REQUIRED)

# Include subdirectory with project sources
add_subdirectory(src)
