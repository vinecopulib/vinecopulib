cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project(vinecopulib VERSION 1.0.0)

enable_testing()

include(cmake/options.cmake               REQUIRED)

include(cmake/findEigen3.cmake            REQUIRED)

include(cmake/findNlopt.cmake             REQUIRED)

include(cmake/findR.cmake                 REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findDoxygen.cmake           REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/paths.cmake                 REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)

include(cmake/gitRevision.cmake           REQUIRED)
