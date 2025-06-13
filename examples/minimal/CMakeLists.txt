cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)

project (Example)

if(NOT DEFINED VINECOPULIB_INCLUDE_DIR)
  message(FATAL_ERROR "VINECOPULIB_INCLUDE_DIR must be specified (e.g., -DVINECOPULIB_INCLUDE_DIR=/path/to/include)")
else()
  message(STATUS "Using VINECOPULIB_INCLUDE_DIR: ${VINECOPULIB_INCLUDE_DIR}")
endif()
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  message(FATAL_ERROR "EIGEN3_INCLUDE_DIR must be specified (e.g., -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3)")
else()
  message(STATUS "Using EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
endif()
if(NOT DEFINED BOOST_INCLUDE_DIR)
  message(FATAL_ERROR "BOOST_INCLUDE_DIR must be specified (e.g., -DBOOST_INCLUDE_DIR=/usr/include)")
else()
  message(STATUS "Using BOOST_INCLUDE_DIR: ${BOOST_INCLUDE_DIR}")
endif()


include_directories(
  "${VINECOPULIB_INCLUDE_DIR}"
)
include_directories(
  SYSTEM
  "${EIGEN3_INCLUDE_DIR}"
  "${BOOST_INCLUDE_DIR}"
)
add_executable(main main.cpp)