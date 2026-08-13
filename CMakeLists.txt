# 3.14 for FetchContent_MakeAvailable; do not lower.
cmake_minimum_required(VERSION 3.14)

project(vinecopulib VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(cmake/options.cmake               REQUIRED)

if(VINECOPULIB_IS_TOP_LEVEL AND BUILD_TESTING)
  enable_testing()
endif()

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findHeaders.cmake           REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)
