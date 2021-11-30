include(cmake/findEigen3.cmake            REQUIRED)
message(STATUS "Eigen3 version: ${EIGEN3_VERSION}")
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)
find_package(wdm                          REQUIRED)

if(TORCH_INSTALL_PREFIX)
  find_package(Torch REQUIRED PATHS ${TORCH_INSTALL_PREFIX})
endif()

set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})
