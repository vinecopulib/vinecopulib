project(Test)

project(Example)

include(CTest)
enable_testing()

#set(gtest_disable_pthreads on) #needed in MinGW
include(../cmake/DownloadProject.cmake)
download_project(
        PROJ                googletest
        GIT_REPOSITORY      https://github.com/google/googletest.git
        GIT_TAG             master
        UPDATE_DISCONNECTED 1
)

include_directories(${test_includes} ${external_includes}) #
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
add_subdirectory(src_test)

set(unit_tests
        test_bicop_parametric
        test_bicop_class)

foreach(app ${unit_tests})
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} gtest gmock vinecopulib ${external_libs} ${test_libs} src_test)
endforeach()

#add_executable(runUnitTests gtest.cpp) #${test_SRCS}
#target_link_libraries(runUnitTests gtest gmock common ${external_libs} ${test_libs} src_test)#
