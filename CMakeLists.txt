cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project(vinecopulib VERSION 0.2.2)

include(cmake/options.cmake               REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findDoxygen.cmake           REQUIRED)

include(cmake/findHeaders.cmake           REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)
