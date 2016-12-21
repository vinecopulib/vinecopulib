cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project(vinecoplib)

include(cmake/options.cmake               REQUIRED)

include(cmake/findGSL.cmake               REQUIRED)

include(cmake/findEigen3.cmake            REQUIRED)

include(cmake/findR.cmake                 REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findTestDependencies.cmake  REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/paths.cmake                 REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)
