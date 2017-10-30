include(cmake/findEigen3.cmake            REQUIRED)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
