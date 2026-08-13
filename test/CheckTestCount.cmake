# Asserts that test_all registers a plausible number of tests. Guards against a
# link configuration in which the self-registering tests are dropped: the binary
# still builds and exits 0, having run nothing.
execute_process(COMMAND "${TEST_ALL}" --gtest_list_tests
                OUTPUT_VARIABLE listing RESULT_VARIABLE status)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "could not list tests: ${TEST_ALL} exited ${status}")
endif()
string(REGEX MATCHALL "\n  [A-Za-z_]" tests "${listing}")
list(LENGTH tests count)
if(count LESS MINIMUM)
    message(FATAL_ERROR
            "only ${count} tests registered, expected at least ${MINIMUM}; the "
            "test object libraries are probably not linked whole (they must be "
            "OBJECT, not STATIC)")
endif()
message(STATUS "${count} tests registered")
