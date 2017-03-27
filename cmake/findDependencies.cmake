#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.63 REQUIRED)

set(external_includes ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})

set(external_libs "-lpthread" ${NLOPT_LIBRARIES} ${Boost_LIBRARIES})
