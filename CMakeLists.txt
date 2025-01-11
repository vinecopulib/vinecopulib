cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)

project(vinecopulib VERSION 0.7.0)

include(CTest)

include(cmake/options.cmake               REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findHeaders.cmake           REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)
