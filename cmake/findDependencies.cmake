#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.55 REQUIRED)

set(external_includes ${GSL_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

set(external_libs "-lpthread" ${GSL_LIBRARIES} ${NLOPT_LIBRARIES} ${Boost_LIBRARIES})
