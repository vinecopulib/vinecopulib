if (WIN32)
    set(gtest_disable_pthreads on) #needed in MinGW
    set(gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" FORCE)
endif ()

include(../cmake/DownloadProject.cmake)
download_project(
        PROJ googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG master
        UPDATE_DISCONNECTED 1
        QUIET 1
)

include_directories(SYSTEM ${external_includes})
include_directories(${CMAKE_BINARY_DIR}/generated/test/)
add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
add_subdirectory(src_test)

foreach (app ${unit_tests})
    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} gtest vinecopulib ${CMAKE_THREAD_LIBS_INIT} ${test_libs} src_test)
endforeach ()
