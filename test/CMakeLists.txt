if (WIN32)
    set(gtest_disable_pthreads on) #needed in MinGW
endif ()

include_directories(SYSTEM ${external_includes})
include_directories(${CMAKE_BINARY_DIR}/generated/test/)
add_subdirectory(src_test)


include(GoogleTest)
foreach (app ${unit_tests})
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} gtest vinecopulib ${CMAKE_THREAD_LIBS_INIT} ${test_libs} src_test)
endforeach ()
if(NOT WIN32)
  gtest_discover_tests(test_all)
else()
  # Fallback: manually register the test executable as a single CTest entry
  add_test(NAME vinecopulib_tests COMMAND test_all)
endif()
# add_test(NAME ${app} COMMAND ${app})