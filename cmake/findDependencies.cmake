#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

include(cmake/findEigen3.cmake            REQUIRED)
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

find_package(NLopt                        QUIET)
if ( NLOPT_INCLUDE_DIRS )
    # this is needed for GCC 6 because we later use these as SYSTEM headers
    # see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=70129
    #     https://issues.apache.org/jira/browse/THRIFT-3828
    #     https://gitlab.kitware.com/cmake/cmake/issues/16291
    get_filename_component(NLOPT_INCLUDE_DIRS ${NLOPT_INCLUDE_DIRS} REALPATH)
else ( NLOPT_INCLUDE_DIRS )
    # for NLopt versions installed with apt-get/brew
    include(cmake/findNLopt.cmake            REQUIRED)
endif()

set(STAN_VERSION "2.16.0")
set(STAN_MATH_VERSION "2.16.0")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/findStan.cmake)

set(external_includes
        ${EIGEN3_INCLUDE_DIR}
        ${NLOPT_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${STAN_INCLUDE_DIRS})

set(external_libs ${CMAKE_THREAD_LIBS_INIT} ${NLOPT_LIBRARIES} ${Boost_LIBRARIES})
