cmake_minimum_required(VERSION 3.2)
enable_testing()

set(CMAKE_CXX_STANDARD 11)

project(vinecopulib)

include(cmake/options.cmake               REQUIRED)

include(cmake/findEigen3.cmake            REQUIRED)

include(cmake/findNlopt.cmake             REQUIRED)

include(cmake/findR.cmake                 REQUIRED)

include(cmake/findDependencies.cmake      REQUIRED)

include(cmake/findTestDependencies.cmake  REQUIRED)

include(cmake/compilerDefOpt.cmake        REQUIRED)

include(cmake/paths.cmake                 REQUIRED)

include(cmake/buildTargets.cmake          REQUIRED)

include(cmake/printInfo.cmake             REQUIRED)

include(cmake/findDoxygen.cmake           REQUIRED)


# generate a header file with git revision id
if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
  add_custom_target(git_revision.h
   git log -1 "--format=format:#define GIT_REVISION \"%H\"%n" HEAD > include/git_revision.h
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} VERBATIM
  )
endif()
