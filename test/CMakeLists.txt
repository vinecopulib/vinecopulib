if(WIN32)
    set(gtest_disable_pthreads on) # needed in MinGW
endif()

include_directories(SYSTEM ${external_includes})
include_directories(${CMAKE_BINARY_DIR}/generated/test/)
add_subdirectory(src_test)

include(GoogleTest)

# test_all is the suite CI and coverage run.
add_executable(test_all test_all.cpp)
target_link_libraries(test_all PRIVATE gtest vinecopulib Threads::Threads
                                       src_test ${test_areas})

# One binary per area for iterating on that area alone: `make test_bicop_kernel`
# compiles only that area, and the objects are shared with test_all, so nothing
# is compiled twice.
foreach(area ${test_areas})
    add_executable(${area}_only EXCLUDE_FROM_ALL test_all.cpp)
    target_link_libraries(${area}_only PRIVATE gtest vinecopulib
                                               Threads::Threads src_test ${area})
endforeach()

# PRE_TEST defers enumeration to test time, so cross-compiled and
# not-yet-built trees still configure.
gtest_discover_tests(test_all DISCOVERY_MODE PRE_TEST)

# A test registers itself from a static initializer that nothing references, so
# a regression to STATIC object libraries links and runs clean while silently
# registering no tests at all. Fail loudly instead.
add_test(NAME test_discovery_is_not_empty
         COMMAND ${CMAKE_COMMAND}
                 -D "TEST_ALL=$<TARGET_FILE:test_all>"
                 -D "MINIMUM=100"
                 -P "${CMAKE_CURRENT_SOURCE_DIR}/CheckTestCount.cmake")
