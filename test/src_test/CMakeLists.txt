file(GLOB_RECURSE src_test_sources *.cpp *.cc *c)
file(GLOB_RECURSE src_test_headers include/*.h include/*.hpp)

include_directories(${external_includes} ${CMAKE_BINARY_DIR}/generated/test/)

add_library(src_test STATIC ${src_test_sources} ${src_test_headers})
target_link_libraries(src_test gtest vinecopulib ${test_libs} ${CMAKE_THREAD_LIBS_INIT})
